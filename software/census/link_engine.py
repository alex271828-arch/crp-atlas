"""Resumable whole-link growth, with rooted intermediate states and exact pruning."""

import os
import pickle
import time
from pathlib import Path

import boundary_filter as boundary
import link_growth as growth
import surface_pruning as surface
import certified_geometry as geometry

VERSION = 116


def normalized_menu(menu):
    return tuple(sorted({boundary.canonical(tuple(word)) for word in menu}))


def initial(seed, menu):
    return dict(
        version=VERSION,
        menu=normalized_menu(menu),
        agenda=[("outer", frozenset(map(tuple, seed)), None, None)],
        seen_outer=set(),
        seen_inner=set(),
        exhausted=False,
        stats=dict(
            t=0.0,
            outer_states=0,
            inner_states=0,
            faces_generated=0,
            local_kills=0,
            boundary_kills=0,
            prune_reasons={},
            pinched_states=0,
            closed_candidates=0,
            canonical_seconds=0.0,
            filter_seconds=0.0,
            generation_seconds=0.0,
            geometry_seconds=0.0,
            geometry_kills=0,
            geometry_boxes=0,
        ),
    )


def save(state, path):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(state, handle, pickle.HIGHEST_PROTOCOL)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def load(path, menu):
    """Load only this program's trusted local checkpoint, never an untrusted pickle."""
    with open(path, "rb") as handle:
        state = pickle.load(handle)
    if state["version"] != VERSION or state["menu"] != normalized_menu(menu):
        raise ValueError("Checkpoint version/menu mismatch; use a fresh cache")
    return state


def expand(item, state, candidate_sink):
    """Compute all children before committing the agenda/seen update.

    Each interval cover encloses every realization of its labelled fragment.
    Thus an isomorphic representative may safely keep either cover. Inner
    deduplication additionally preserves the selected pivot.
    """
    kind, faces, pivot, region = item
    stats, menu = state["stats"], state["menu"]
    started = time.monotonic()
    info = growth.analyze(faces, menu)
    if info is None:
        stats["local_kills"] += 1
        stats["filter_seconds"] += time.monotonic() - started
        return [], None, None
    verdict = surface.propagate(faces, info, menu)
    stats["filter_seconds"] += time.monotonic() - started
    if verdict["excluded"]:
        stats["boundary_kills"] += 1
        reasons = stats["prune_reasons"]
        reason = verdict["reason"]
        reasons[reason] = reasons.get(reason, 0) + 1
        return [], None, None
    faces, info = verdict["faces"], verdict["info"]
    started = time.monotonic()
    code = growth.canonical_code(faces, pivot if kind == "inner" else None)
    stats["canonical_seconds"] += time.monotonic() - started
    seen = state["seen_inner"] if kind == "inner" else state["seen_outer"]
    if code in seen:
        return [], None, None
    stats["pinched_states"] += bool(info["pinched"])
    if kind == "inner":
        stats["inner_states"] += 1
        if pivot not in info["opened"]:
            children = [("outer", faces, None, region)]
        else:
            started = time.monotonic()
            children = [
                ("inner", faces | {face}, pivot, region)
                for face in growth.next_faces(faces, pivot, menu, info, stats)
            ]
            stats["generation_seconds"] += time.monotonic() - started
    else:
        stats["outer_states"] += 1
        closed = len(info["vertices"]) - len(info["opened"])
        # Only proof effort changes; expiration retains all unrefuted boxes.
        budget = (
            0.05 if closed < 4 else 1.0 if closed < 8 else 3.0 if closed < 12 else 10.0
        )
        started = time.monotonic()
        excluded, region, boxes = geometry.solve(faces, info, region, time_cap=budget)
        stats["geometry_seconds"] += time.monotonic() - started
        stats["geometry_boxes"] += boxes
        children = []
        if excluded:
            stats["geometry_kills"] += 1
        elif info["opened"]:
            children = [("inner", faces, growth.choose_pivot(info), region)]
        elif len(info["vertices"]) - len(info["directed"]) // 2 + len(faces) == 2:
            candidate_sink(faces, code)  # durable before this item can be forgotten
            stats["closed_candidates"] += 1
    return children, seen, code


def run(state, seconds, candidate_sink, checkpoint=None, stop=None, max_steps=None):
    """Pause between complete expansions; never checkpoint half of a sibling set."""
    started = accounted = last_saved = time.monotonic()
    steps = 0
    while state["agenda"] and time.monotonic() - started < seconds:
        if (stop and stop()) or (max_steps is not None and steps >= max_steps):
            break
        # Leave the item on the agenda until expansion and journalling succeed.
        children, seen, code = expand(state["agenda"][-1], state, candidate_sink)
        state["agenda"].pop()
        state["agenda"].extend(children)
        if seen is not None:
            seen.add(code)
        steps += 1
        now = time.monotonic()
        if checkpoint and now - last_saved >= 20:
            state["stats"]["t"] += now - accounted
            accounted = last_saved = now
            save(state, checkpoint)
    state["stats"]["t"] += time.monotonic() - accounted
    state["exhausted"] = not state["agenda"]
    if checkpoint:
        save(state, checkpoint)
    return state
