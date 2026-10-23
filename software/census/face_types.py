"""Exact face angles and curvature units for the finite face alphabet."""

from functools import lru_cache

FACE_SIZES = (3, 4, 5, 6, 8, 10)
ANGLE_UNIT_DEGREES = 3
FULL_TURN = 120
HALF_TURN = 60
TOTAL_CURVATURE = 240
ANGLE = {n: 60 * (n - 2) // n for n in FACE_SIZES}
ANGLE_DEGREES = {n: ANGLE_UNIT_DEGREES * a for n, a in ANGLE.items()}
assert all(60 * (n - 2) % n == 0 for n in FACE_SIZES)


@lru_cache(None)
def defect(word):
    """Exact vertex curvature in three-degree units, for a cyclic face word."""
    return FULL_TURN - sum(ANGLE[n] for n in word)
