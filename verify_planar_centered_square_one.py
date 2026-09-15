#!/usr/bin/env python3
"""Verify the planar W_4(1) four-copy inequality."""

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import sqrt

import mpmath as mp


mp.iv.dps = 260
mp.mp.dps = 60


@dataclass(frozen=True)
class Quadratic:
    rational: F
    radical: F

    def __add__(self, other):
        return Quadratic(self.rational + other.rational,
                         self.radical + other.radical)

    def __neg__(self):
        return Quadratic(-self.rational, -self.radical)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Quadratic(
            self.rational * other.rational + 2 * self.radical * other.radical,
            self.rational * other.radical + self.radical * other.rational,
        )


ZERO = Quadratic(F(0), F(0))
ONE = Quadratic(F(1), F(0))
RADIUS = Quadratic(F(0), F(1, 2))
CENTERS = (
    (ZERO, ZERO), (ONE, ZERO), (ONE, ONE), (ZERO, ONE),
)
POINT_OFFSETS = (
    (ZERO, ZERO), (RADIUS, ZERO), (ZERO, RADIUS),
    (-RADIUS, ZERO), (ZERO, -RADIUS),
)
POINTS = tuple(
    (x + dx, y + dy)
    for x, y in CENTERS
    for dx, dy in POINT_OFFSETS
)
PAIRS = tuple(combinations(range(len(POINTS)), 2))
DISTANCES = {
    pair: sum((
        (POINTS[pair[0]][coordinate] - POINTS[pair[1]][coordinate])
        * (POINTS[pair[0]][coordinate] - POINTS[pair[1]][coordinate])
        for coordinate in (0, 1)
    ), ZERO)
    for pair in PAIRS
}
EDGES = tuple(pair for pair in PAIRS if DISTANCES[pair] == ONE)
BLOCKS = tuple(tuple(5 * block + i for i in range(5)) for block in range(4))

LOCAL_SHIFT = F(499999, 500000)
CONSTANT = F(-23994231, 1000000) - LOCAL_SHIFT
LINEAR_ORBITS = (
    ((0, 5, 10, 15), F("3.167237")),
    ((1, 2, 7, 8, 13, 14, 16, 19), F("5.837166")),
    ((3, 4, 6, 9, 11, 12, 17, 18), F(10)),
)
LINEAR_BY_INDEX = {
    index: value for orbit, value in LINEAR_ORBITS for index in orbit
}
LINEAR = tuple(LINEAR_BY_INDEX.get(index, F(0)) for index in range(20))

QUADRATIC_ORBITS = (
    ((0, 1, 86, 87, 147, 148, 180, 183), F("3.129119")),
    ((2, 3, 85, 88, 145, 146, 181, 182), F(-10)),
    ((4, 14, 89, 149), F(10)),
    ((5, 16, 55, 80, 91, 135, 152, 157), F("-2.904617")),
    ((6, 15, 32, 39, 92, 125, 153, 175), F("0.355813")),
    ((7, 18, 22, 49, 93, 114, 150, 170), F("-2.904617")),
    ((8, 17, 65, 70, 90, 102, 151, 164), F("-0.618433")),
    ((9, 94), F(10)),
    ((10, 11, 60, 75, 96, 97, 107, 140), F("-0.487123")),
    ((12, 13, 27, 44, 95, 98, 119, 130), F("9.381567")),
    ((19, 112, 169, 186), F(10)),
    ((20, 38, 100, 113, 155, 163, 185, 188), F("-7.095383")),
    ((21, 37, 99, 124, 156, 162, 184, 189), F(10)),
    ((23, 51, 58, 84, 116, 139, 158, 173), F(10)),
    ((24, 36, 42, 50, 117, 129, 174, 176), F(0)),
    ((25, 53, 118, 171), F("-8.479162")),
    ((26, 52, 69, 73, 106, 115, 165, 172), F("-0.487123")),
    ((28, 46, 63, 79, 108, 121, 133, 144), F("-0.487123")),
    ((29, 45, 64, 78, 111, 122, 132, 141), F("-0.487123")),
    ((30, 48, 123, 131), F("4.294145")),
    ((31, 47, 120, 134), F("4.162834")),
    ((33, 41, 128, 179), F(10)),
    ((34, 40, 57, 81, 127, 138, 161, 178), F(0)),
    ((35, 43, 66, 72, 105, 126, 168, 177), F("2.286183")),
    ((54, 101, 154, 187), F(10)),
    ((56, 82, 137, 160), F("0.355813")),
    ((59, 67, 71, 83, 104, 136, 159, 167), F("0.487123")),
    ((61, 77, 110, 142), F(0)),
    ((62, 76, 109, 143), F("0.131310")),
    ((68, 74, 103, 166), F(10)),
)
QUADRATIC_BY_INDEX = {
    index: value for orbit, value in QUADRATIC_ORBITS for index in orbit
}
assert len(QUADRATIC_BY_INDEX) == len(PAIRS)
QUADRATIC = tuple(QUADRATIC_BY_INDEX[index] for index in range(len(PAIRS)))

LOCAL_FLOOR = F(0)
LAMBDA = F(46)
MU = F(600)
TAIL_START = F(125)
SERIES_DEGREE_FLOOR = 180


def allowed_masks():
    return [
        mask for mask in range(1 << len(POINTS))
        if not any(
            mask >> i & 1 and mask >> j & 1
            for i, j in EDGES
        )
        and all(any(mask >> i & 1 for i in block) for block in BLOCKS)
    ]


def local_value(mask):
    return (
        CONSTANT
        + sum(LINEAR[i] for i in range(len(POINTS)) if mask >> i & 1)
        + sum(
            coefficient
            for coefficient, (i, j) in zip(QUADRATIC, PAIRS)
            if mask >> i & 1 and mask >> j & 1
        )
    )


def check_local():
    masks = allowed_masks()
    values = [local_value(mask) for mask in masks]
    assert len(EDGES) == 36
    assert len(masks) == 2338
    assert min(values) == LOCAL_FLOOR
    print(f"local: {len(masks)} masks; minimum {min(values)}")


def check_density():
    linear = sum(LINEAR)
    quadratic = sum(QUADRATIC)
    assert CONSTANT == -F(24994229, 1000000)
    assert linear == F(34841569, 250000)
    assert quadratic == F(37900733, 125000)
    assert quadratic - LAMBDA - MU == -F(42849267, 125000)
    coefficient = quadratic - LAMBDA - MU
    linear += LAMBDA
    candidates = [F(1, 5), F(1, 4)]
    if coefficient:
        vertex = -linear / (2 * coefficient)
        if F(1, 5) < vertex < F(1, 4):
            candidates.append(vertex)
    values = [
        CONSTANT + linear * p + coefficient * p * p
        for p in candidates
    ]
    maximum = max(values)
    assert maximum == -F(154587, 2000000)
    assert maximum < LOCAL_FLOOR
    print(f"coefficients: A={linear - LAMBDA}, S={quadratic}, S-646={coefficient}")
    print(f"density: maximum {maximum} < {LOCAL_FLOOR}")


def fraction_interval(value):
    return mp.iv.mpf(value.numerator) / value.denominator


def algebraic_interval(value):
    rational = fraction_interval(value.rational)
    radical = fraction_interval(value.radical)
    return rational + radical * mp.iv.sqrt(mp.iv.mpf(2))


def series_degree(z):
    upper = float(mp.mpf(z.b))
    return max(SERIES_DEGREE_FLOOR, int(2.2 * sqrt(max(0.0, upper)) + 30))


def j0_interval(q, midpoint):
    q_iv = algebraic_interval(q)
    t_iv = fraction_interval(midpoint)
    z = q_iv * t_iv * t_iv / 4
    degree = series_degree(z)
    term = mp.iv.mpf(1)
    partial = term
    for index in range(1, degree + 1):
        term *= z / (index * index)
        partial = partial - term if index & 1 else partial + term
    target = mp.iv.mpf("1e-80")
    while True:
        ratio = z / ((degree + 2) ** 2)
        if ratio < 1:
            tail = term * z / ((degree + 1) ** 2) / (1 - ratio)
            if tail < target:
                return partial + mp.iv.mpf([-1, 1]) * tail
        for index in range(degree + 1, degree + 21):
            term *= z / (index * index)
            partial = partial - term if index & 1 else partial + term
        degree += 20


def j1_interval(q, midpoint):
    q_iv = algebraic_interval(q)
    t_iv = fraction_interval(midpoint)
    z = q_iv * t_iv * t_iv / 4
    degree = series_degree(z)
    term = mp.iv.mpf(1)
    partial = term
    for index in range(1, degree + 1):
        term *= z / (index * (index + 1))
        partial = partial - term if index & 1 else partial + term
    scale = mp.iv.sqrt(q_iv) * t_iv / 2
    target = mp.iv.mpf("1e-80")
    while True:
        ratio = z / ((degree + 2) * (degree + 3))
        if ratio < 1:
            tail = (term * z / ((degree + 1) * (degree + 2))
                    / (1 - ratio)) * scale
            if tail < target:
                return partial * scale + mp.iv.mpf([-1, 1]) * tail
        for index in range(degree + 1, degree + 21):
            term *= z / (index * (index + 1))
            partial = partial - term if index & 1 else partial + term
        degree += 20


def kernel_terms():
    terms = {}
    for distance, coefficient in zip(
        (DISTANCES[pair] for pair in PAIRS), QUADRATIC
    ):
        terms[distance] = terms.get(distance, F(0)) + coefficient
    terms[ONE] -= MU
    return tuple((distance, coefficient) for distance, coefficient in terms.items()
                 if coefficient)


def cell_upper(terms, lower_t, upper_t):
    midpoint = (lower_t + upper_t) / 2
    half_width = fraction_interval((upper_t - lower_t) / 2)
    value = mp.iv.mpf(0)
    derivative = mp.iv.mpf(0)
    second_derivative = mp.iv.mpf(0)
    for distance, coefficient in terms:
        coefficient_iv = fraction_interval(coefficient)
        distance_iv = algebraic_interval(distance)
        value += coefficient_iv * j0_interval(distance, midpoint)
        derivative -= (coefficient_iv * mp.iv.sqrt(distance_iv)
                       * j1_interval(distance, midpoint))
        second_derivative += abs(coefficient_iv) * distance_iv
    return (value + abs(derivative) * half_width
            + second_derivative * half_width ** 2 / 2)


def check_cells(terms, lower, upper, denominator):
    maximum = float("-inf")
    maximum_cell = None
    threshold = fraction_interval(LAMBDA)
    first = int(lower * denominator)
    last = int(upper * denominator)
    for index in range(first, last):
        lower_t = F(index, denominator)
        upper_t = F(index + 1, denominator)
        value = cell_upper(terms, lower_t, upper_t)
        assert value < threshold
        display_upper = float(mp.mpf(value.b))
        if display_upper > maximum:
            maximum = display_upper
            maximum_cell = (lower_t, upper_t)
    return maximum, maximum_cell


def check_kernel():
    terms = kernel_terms()
    compact, cell = check_cells(terms, F(0), TAIL_START, 20)
    fourth_root_sum = sum(
        fraction_interval(abs(coefficient))
        / algebraic_interval(distance) ** mp.iv.mpf("0.25")
        for distance, coefficient in terms
    )
    assert fourth_root_sum < fraction_interval(F(542))
    tail = mp.iv.sqrt(fraction_interval(F(2, 3))
                      / fraction_interval(TAIL_START)) * fourth_root_sum
    assert tail < fraction_interval(F(40))
    print(f"kernel: maximum cell upper {compact:.12f} at {cell}")
    print(f"kernel: {len(terms)} terms; every interval of width 1/20 has H<46")
    print(f"kernel: tail < {mp.nstr(mp.mpf(tail.b), 12)} < 40")


def main():
    check_local()
    check_density()
    check_kernel()
    print(r"E^2 ->_m (ell_2, W_4(1)) passed")


if __name__ == "__main__":
    main()
