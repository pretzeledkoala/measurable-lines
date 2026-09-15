#!/usr/bin/env python3
"""Verify the finite local data for ``(ell_2, ell_7)``."""

import collections
from fractions import Fraction
from itertools import combinations


NV = 18
LOCAL_SHIFT = -1
PAPER_LOCAL_BOUND = 0
VERTEX_COEFF = (
    0, -54, -54, -16, -2, 15, -46, -1, 26, 18, -16, -52,
    1, -18, -23, -44, -53, 0,
)
EDGE_COEFF = (
    -100, 0, 0, 0, 0, 0, 0, 0, -100, 0, 0, 0,
    0, 0, 0, 0, 0, -100, -2, -11, -12, 14, -16, 11,
    -100, -100, 40, 14, 31, 30, 57, 26, 0, -100, -54, -100,
    -58, -22, -43, 77, -100, -100, 99, 78, 79, 100, 59, 0,
    -100, -95, -100, -60, -55, 63, 85, -100, -100, 93, 67, 68,
    26, 0, -100, -87, -92, -58, 63, 57, 100, -100, -100, 84,
    56, 28, 0, -100, -63, -58, 52, 61, 65, 70, -100, -100,
    85, 7, 0, -100, -42, 67, 96, 100, 73, 100, -89, 100,
    39, 0, 100, -85, 93, 90, 58, 60, 88, 100, 100, 0,
    -90, -80, 63, 47, 60, 60, 75, 100, 100, 100, -38, -52,
    -52, -49, -43, 5, 0, 100, -57, -77, -46, -13, -12, 0,
    100, -77, -100, -57, 20, 0, 100, -95, -100, -12, 0, 100,
    -60, -10, 0, 100, 1, 0, 100, 0, 100,
)
KERNEL_Q = (1, 3, 4, 7, 9, 13, 16, 21, 25, 31, 36, 43, 49, 57)
KERNEL_COEFF = (
    -489, 1035, -669, 600, -700, 548, -355,
    441, -105, 339, -114, -139, 16, -90,
)
SQRT_UP_NUM = (
    1000000, 1732051, 2000000, 2645752, 3000000, 3605552,
    4000000, 4582576, 5000000, 5567765, 6000000, 6557439,
    7000000, 7549835,
)
FOURTH_LB_NUM = (
    1, 13, 7, 8, 17, 189, 2, 21, 11, 47, 12, 5, 13, 27,
)
FOURTH_LB_DEN = (
    1, 10, 5, 5, 10, 100, 1, 10, 5, 20, 5, 2, 5, 10,
)


def qform(a, b):
    ia, ja = divmod(a, 9)
    ib, jb = divmod(b, 9)
    x, y = ia - ib, ja - jb
    return x * x + x * y + y * y


EDGES = tuple(
    (a, b, qform(a, b)) for a, b in combinations(range(NV), 2)
)


def verify_data():
    assert len(VERTEX_COEFF) == NV
    assert len(EDGES) == len(EDGE_COEFF) == 153

    aggregates = collections.Counter()
    for (a, b, distance), coefficient in zip(EDGES, EDGE_COEFF):
        assert distance == qform(a, b)
        aggregates[distance] += coefficient

    assert sum(VERTEX_COEFF) == -319
    assert sum(abs(c) for c in VERTEX_COEFF) == 439
    assert sum(abs(c) for c in EDGE_COEFF) == 8562
    red_linear = [Fraction(500, 9) + c for c in VERTEX_COEFF]
    assert sum(red_linear) == 681
    assert sum(abs(c) for c in red_linear) == 681
    assert 18 * Fraction(500, 9) == 1000

    for distance, coefficient in zip(KERNEL_Q, KERNEL_COEFF):
        assert aggregates[distance] == coefficient
    assert aggregates[64] == aggregates[73] == 0
    assert -319 + sum(KERNEL_COEFF) == -1

    derivative_sum = sum(
        Fraction(abs(c) * n, 10**6)
        for c, n in zip(KERNEL_COEFF, SQRT_UP_NUM)
    )
    assert all(n * n >= q * 10**12 for q, n in zip(KERNEL_Q, SQRT_UP_NUM))
    assert Fraction(2, 3) * derivative_sum == Fraction(35046648006, 3000000)

    tail_sum = sum(
        Fraction(abs(c) * d, p)
        for c, p, d in zip(KERNEL_COEFF, FOURTH_LB_NUM, FOURTH_LB_DEN)
    )
    assert all(
        p**4 <= q * d**4
        for q, p, d in zip(KERNEL_Q, FOURTH_LB_NUM, FOURTH_LB_DEN)
    )
    assert tail_sum == Fraction(384576946829, 107972865) < 3600

    multiplicities = collections.Counter(distance for _, _, distance in EDGES)
    assert multiplicities == {
        1: 33,
        3: 15,
        4: 14,
        7: 13,
        9: 12,
        13: 11,
        16: 10,
        21: 9,
        25: 8,
        31: 7,
        36: 6,
        43: 5,
        49: 4,
        57: 3,
        64: 2,
        73: 1,
    }

def allowed(mask):
    if any(
        distance == 1
        and (mask >> a) & 1
        and (mask >> b) & 1
        for a, b, distance in EDGES
    ):
        return False

    return all(
        any((mask >> (9 * column + j)) & 1 for j in range(start, start + 7))
        for column in range(2)
        for start in range(3)
    )


def functional(mask):
    value = LOCAL_SHIFT + sum(
        coefficient
        for vertex, coefficient in enumerate(VERTEX_COEFF)
        if (mask >> vertex) & 1
    )
    return value + sum(
        coefficient
        for (a, b, _), coefficient in zip(EDGES, EDGE_COEFF)
        if (mask >> a) & 1 and (mask >> b) & 1
    )


def verify_local():
    accepted = 0
    minimum = None
    minimizers = 0
    by_size = {}

    for mask in range(1 << NV):
        if not allowed(mask):
            continue
        accepted += 1
        size = mask.bit_count()
        value = functional(mask)
        assert value >= PAPER_LOCAL_BOUND, (mask, value)
        count, current = by_size.get(size, (0, None))
        by_size[size] = (
            count + 1,
            value if current is None else min(current, value),
        )
        if minimum is None or value < minimum:
            minimum = value
            minimizers = 1
        elif value == minimum:
            minimizers += 1

    assert accepted == 751
    assert minimum == PAPER_LOCAL_BOUND
    assert minimizers == 36
    assert by_size == {
        2: (16, 1),
        3: (150, 0),
        4: (337, 0),
        5: (220, 0),
        6: (28, 0),
    }

def main():
    verify_data()
    verify_local()


if __name__ == "__main__":
    main()
