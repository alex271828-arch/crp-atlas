"""Pass partition: octagon/decagon links first, then the unchanged v1.11 order."""

from itertools import product
from boundary_filter import canonical
from face_types import FACE_SIZES, ANGLE_DEGREES as ANG


def fan(link):
    """Central fan realizing a cyclic link of face sizes. Apex = 0; spokes and interior
    vertices are labelled 1..M sequentially around the rim. Returns (faces, rim_set).
    A face of size n spans spoke s0, its (n-3) fresh interior vertices, then spoke s1.
    """
    d = len(link)
    faces = []
    label = 1
    for i in range(d):
        n = link[i]
        s0 = label
        interiors = [label + 1 + t for t in range(n - 3)]
        label_next = label + (n - 2)
        s1 = 1 if i == d - 1 else label_next  # last face wraps to spoke 1
        faces.append(tuple([0, s0] + interiors + [s1]))
        label = label_next
    M = sum(n - 2 for n in link)
    return faces, set(range(1, M + 1))


def name(link):
    """Readable unique name: 'a^k' for uniform links, else dotted canonical (e.g. 3.4.3.4)."""
    if len(set(link)) == 1:
        return "%d^%d" % (link[0], len(link))
    return ".".join(str(x) for x in link)


def enumerate_links(menu=FACE_SIZES):
    return sorted(
        {
            canonical(w)
            for d in range(3, 6)
            for w in product(menu, repeat=d)
            if sum(ANG[n] for n in w) < 360
        },
        key=lambda s: (len(s), s.count(6), s.count(5), s.count(4), s),
    )


def baseline_order(words):
    return sorted(
        words,
        key=lambda s: (
            int(s.count(3) >= 3),
            -(360 - sum(ANG[n] for n in s)),
            len(s),
            (len(s), s.count(6), s.count(5), s.count(4), s),
        ),
    )


ALL_LINKS = enumerate_links()
OLD_LINKS = baseline_order([w for w in ALL_LINKS if max(w) <= 6])
NEW_LINKS = baseline_order([w for w in ALL_LINKS if max(w) > 6])
SIGS = NEW_LINKS + OLD_LINKS
ORDER = [name(w) for w in SIGS]
SEEDS = {}
for w, label in zip(SIGS, ORDER):
    f, rim = fan(w)
    SEEDS[label] = dict(central=f, rim=rim, link=w)
BASELINE_OFFSET = len(NEW_LINKS)
BASELINE_PASS = {k: k + BASELINE_OFFSET for k in range(len(OLD_LINKS))}
