"""Complete face-complex automorphisms, using only polygon combinatorics."""

import argparse
import json
from pathlib import Path
from metric import orient_polygons
from proof import require


def compute(faces):
    """Enumerate automorphisms by mapping one face and propagating across edges.

    A connected oriented polygonal sphere is determined by face-edge adjacency.
    An automorphism is determined by the image of a single directed face edge
    and whether it preserves or reverses orientation. Every possible image is
    tried; each induces at most one bijection, checked on all face cycles.
    This is a numerical optimization only, never a certification assumption.
    """
    faces = orient_polygons(faces)
    vertices = sorted({v for f in faces for v in f})
    require(vertices == list(range(len(vertices))), "Noncontiguous vertex labels")
    across = {}
    for i, face in enumerate(faces):
        for a, b in zip(face, face[1:] + face[:1]):
            across[a, b] = i
    permutations = []
    source = faces[0]
    for target in faces:
        if len(target) != len(source):
            continue
        for direction in (1, -1):
            for offset in range(len(target)):
                mapping, inverse, visited = {}, {}, set()
                pending = [
                    (
                        0,
                        [
                            target[(offset + direction * k) % len(target)]
                            for k in range(len(target))
                        ],
                    )
                ]
                valid = True
                while pending and valid:
                    index, image = pending.pop()
                    face = faces[index]
                    for a, b in zip(face, image):
                        if (a in mapping and mapping[a] != b) or (b in inverse and inverse[b] != a):
                            valid = False
                            break
                        mapping[a], inverse[b] = b, a
                    if not valid or index in visited:
                        continue
                    visited.add(index)
                    for a, b in zip(face, face[1:] + face[:1]):
                        neighbor = across[b, a]
                        if neighbor in visited:
                            continue
                        u, v = mapping[b], mapping[a]
                        next_index = across[u, v] if direction == 1 else across[v, u]
                        next_face = faces[next_index]
                        original = faces[neighbor]
                        if len(next_face) != len(original):
                            valid = False
                            break
                        shift = next_face.index(u) - direction * original.index(b)
                        pending.append(
                            (
                                neighbor,
                                [
                                    next_face[(shift + direction * k) % len(next_face)]
                                    for k in range(len(next_face))
                                ],
                            )
                        )
                if valid and len(visited) == len(faces) and len(mapping) == len(vertices):
                    permutations.append([mapping[v] for v in vertices])
    permutations.sort()
    require(bool(permutations), "Identity automorphism was lost")
    orbits = sorted({tuple(sorted({p[v] for p in permutations})) for v in vertices})
    return dict(
        schema="crp-face-symmetry/1",
        order=len(permutations),
        vertex_orbits=[list(o) for o in orbits],
        permutations=permutations,
        input_faces=faces,
    )


def write_metric(path, faces, orbits):
    lines = ["1", f"{1+max(v for f in faces for v in f)} {len(faces)}"]
    lines.extend(" ".join(map(str, [len(f), *f])) for f in faces)
    lines.append(str(len(orbits)))
    lines.extend(" ".join(map(str, [len(o), *o])) for o in orbits)
    with Path(path).open("x") as stream:
        stream.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("metric", type=Path)
    args = parser.parse_args()
    faces = json.loads(args.input.read_text())["faces"]
    result = compute(faces)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    write_metric(args.metric, faces, result["vertex_orbits"])


if __name__ == "__main__":
    main()
