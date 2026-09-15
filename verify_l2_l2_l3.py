#!/usr/bin/env python3
"""Verify the finite and kernel bounds for ``(ell_2, ell_2, ell_3)``."""

from fractions import Fraction
from itertools import combinations, product


# Colors are encoded by G = 0, R = 1, and B = 2.
G, R, B = range(3)
PAPER_LOCAL_BOUND = 0
POINTS = tuple(product(range(3), repeat=2))
INDEX = {point: index for index, point in enumerate(POINTS)}
DIRECTIONS = ((1, 0), (0, 1), (1, -1))


def q(p, r):
    a, b = p[0] - r[0], p[1] - r[1]
    return a * a + a * b + b * b


EDGES = tuple(
    (p, r) for p, r in combinations(POINTS, 2) if q(p, r) == 1
)
TRIPLES = tuple(
    (p, (p[0] + a, p[1] + b), (p[0] + 2 * a, p[1] + 2 * b))
    for p in POINTS
    for a, b in DIRECTIONS
    if all(
        (p[0] + k * a, p[1] + k * b) in INDEX for k in (1, 2)
    )
)

PATCH_SYMMETRIES = (
    lambda p: p,
    lambda p: (p[1], p[0]),
    lambda p: (2 - p[0], 2 - p[1]),
    lambda p: (2 - p[1], 2 - p[0]),
)


def vertex_orbit(point):
    return sorted({symmetry(point) for symmetry in PATCH_SYMMETRIES})


def pair_orbit(pair):
    p, r = pair
    return sorted({
        tuple(sorted((symmetry(p), symmetry(r))))
        for symmetry in PATCH_SYMMETRIES
    })


VERTEX_DATA = (
    ((0, 0), -252),
    ((0, 1), -1004),
    ((0, 2), 998),
    ((1, 1), -345),
)
PAIR_DATA = (
    ((0, 0), (0, 1), -998, 2),
    ((0, 0), (0, 2), -1000, 0),
    ((0, 0), (1, 1), 500, 150),
    ((0, 0), (1, 2), 1000, 300),
    ((0, 0), (2, 2), 375, 200),
    ((0, 1), (0, 2), -446, -894),
    ((0, 1), (1, 0), 2, 502),
    ((0, 1), (1, 1), 444, 948),
    ((0, 1), (1, 2), 662, -12),
    ((0, 1), (2, 0), 448, 0),
    ((0, 1), (2, 1), -264, -12),
    ((0, 2), (1, 1), -52, -1000),
    ((0, 2), (2, 0), 112, 0),
)
VERTEX_ORBITS = tuple(vertex_orbit(point) for point, _ in VERTEX_DATA)
PAIR_ORBITS = tuple(pair_orbit((p, r)) for p, r, _, _ in PAIR_DATA)


def green(coloring, point):
    return int(coloring[INDEX[point]] == G)


def signed(coloring, point):
    color = coloring[INDEX[point]]
    return 0 if color == G else 1 if color == R else -1


def avoids_forbidden_patterns(coloring):
    no_red_blue_pair = all(
        not (
            coloring[INDEX[p]] == coloring[INDEX[r]]
            and coloring[INDEX[p]] in (R, B)
        )
        for p, r in EDGES
    )
    no_green_triple = all(
        not all(coloring[INDEX[point]] == G for point in triple)
        for triple in TRIPLES
    )
    return no_red_blue_pair and no_green_triple


def four_times_local_function(coloring):
    value = -4 * 196
    for (_, coefficient), orbit in zip(VERTEX_DATA, VERTEX_ORBITS):
        value += coefficient * (4 // len(orbit)) * sum(
            green(coloring, point) for point in orbit
        )
    for (_, _, green_coefficient, signed_coefficient), orbit in zip(
        PAIR_DATA, PAIR_ORBITS
    ):
        scale = 4 // len(orbit)
        value += green_coefficient * scale * sum(
            green(coloring, p) * green(coloring, r) for p, r in orbit
        )
        value += signed_coefficient * scale * sum(
            signed(coloring, p) * signed(coloring, r) for p, r in orbit
        )
    return value


GREEN_KERNEL_TERMS = (
    (-1050, 1),
    (1610, 3),
    (-1152, 4),
    (1000, 7),
    (375, 12),
)
SIGNED_KERNEL_TERMS = (
    (-442, 1),
    (138, 3),
    (-12, 4),
    (300, 7),
    (200, 12),
)
DIGITS = 90
UNIT = 10**DIGITS
TOLERANCE = 10 ** (DIGITS - 35)
N = 400
NODE_COUNT = 10000 + 1
Y_DENOMINATOR = 4 * N * N
SQRT_UPPER = {
    1: Fraction(1),
    3: Fraction(7, 4),
    4: Fraction(2),
    7: Fraction(8, 3),
    12: Fraction(7, 2),
}
ROOT_LOWER = {
    1: Fraction(1),
    3: Fraction(13, 10),
    4: Fraction(7, 5),
    7: Fraction(8, 5),
    12: Fraction(9, 5),
}


def ceiling_division(a, b):
    return -((-a) // b)


def j0_interval(p):
    """Return a directed-rational enclosure for J_0(sqrt(p/640000))."""
    term_lower = term_upper = UNIT
    sum_lower = sum_upper = UNIT
    n = 0
    while True:
        n += 1
        denominator = Y_DENOMINATOR * n * n
        next_lower = term_lower * p // denominator
        next_upper = ceiling_division(term_upper * p, denominator)

        decreasing = p <= Y_DENOMINATOR * (n + 1) ** 2
        if decreasing and next_upper < TOLERANCE:
            if n % 2 == 0:
                sum_upper += next_upper
            else:
                sum_lower -= next_upper
            return sum_lower, sum_upper

        if n % 2 == 0:
            sum_lower += next_lower
            sum_upper += next_upper
        else:
            sum_lower -= next_upper
            sum_upper -= next_lower
        term_lower, term_upper = next_lower, next_upper


def kernel_interval(j, constant, terms):
    lower = upper = constant * UNIT
    for coefficient, distance in terms:
        bessel_lower, bessel_upper = j0_interval(distance * j * j)
        if coefficient >= 0:
            lower += coefficient * bessel_lower
            upper += coefficient * bessel_upper
        else:
            lower += coefficient * bessel_upper
            upper += coefficient * bessel_lower
    return lower, upper


def verify_local():
    valid_colorings = [
        coloring
        for coloring in product((G, R, B), repeat=9)
        if avoids_forbidden_patterns(coloring)
    ]
    assert len(EDGES) == 16
    assert len(TRIPLES) == 7
    assert len(valid_colorings) == 304
    values = [four_times_local_function(coloring) for coloring in valid_colorings]
    assert all(value >= PAPER_LOCAL_BOUND for value in values)
    assert min(values) == PAPER_LOCAL_BOUND


def verify_kernel_nodes():
    for j in range(NODE_COUNT):
        _, green_upper = kernel_interval(j, -603, GREEN_KERNEL_TERMS)
        _, signed_upper = kernel_interval(j, 0, SIGNED_KERNEL_TERMS)
        assert green_upper <= 180 * UNIT
        assert signed_upper < Fraction(46019, 250) * UNIT


def verify_rational_bounds():
    assert Y_DENOMINATOR == 640000
    assert Fraction(10000, N) == 25
    assert all(
        Fraction(distance) <= upper**2
        for distance, upper in SQRT_UPPER.items()
    )
    assert all(
        Fraction(distance) > lower**4
        for distance, lower in ROOT_LOWER.items()
        if distance != 1
    )

    linear_error = sum(abs(coefficient) for _, coefficient in VERTEX_DATA)
    quadratic_error = sum(
        abs(green_coefficient) + abs(signed_coefficient)
        for _, _, green_coefficient, signed_coefficient in PAIR_DATA
    )
    assert linear_error == 2599
    assert quadratic_error == 10323
    assert 4 * linear_error + 12 * quadratic_error == 134272

    green_derivative_bound = sum(
        abs(coefficient) * SQRT_UPPER[distance]
        for coefficient, distance in GREEN_KERNEL_TERMS
    )
    signed_derivative_bound = sum(
        abs(coefficient) * SQRT_UPPER[distance]
        for coefficient, distance in SIGNED_KERNEL_TERMS
    )
    assert Fraction(180) + green_derivative_bound / 800 < 193
    assert Fraction(46019, 250) + signed_derivative_bound / 800 < 193

    assert Fraction(2, 3 * 25) * (
        1 + Fraction(1, 8 * 25) + Fraction(9, 128 * 25 * 25)
    ) ** 2 < Fraction(1, 36)
    assert absolute_tail_bound(GREEN_KERNEL_TERMS, -603) < 55
    assert absolute_tail_bound(SIGNED_KERNEL_TERMS) < 143


def absolute_tail_bound(terms, constant=0):
    return constant + sum(
        Fraction(abs(coefficient), 6) / ROOT_LOWER[distance]
        for coefficient, distance in terms
    )


def main():
    verify_local()
    verify_kernel_nodes()
    verify_rational_bounds()


if __name__ == "__main__":
    main()
