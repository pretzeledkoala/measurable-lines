#!/usr/bin/env python3
"""Verify the square--tetrahedron kernel and analytic bounds."""

from fractions import Fraction
from mpmath import iv


iv.dps = 80
NODE_DENOMINATOR = 1000
NODE_COUNT = 40001
ROOT3 = iv.sqrt(3)
DISTANCE_COEFFICIENTS = (
    (2, -1, -68300),
    (1, 0, -300708),
    (3, -1, 56588),
    (2, 0, -142872),
    (4, -1, -14632),
    (3, 0, 18334),
    (5, -1, 10996),
    (2, 1, 628),
    (4, 0, -5700),
    (3, 1, -404),
    (4, 1, 204),
    (6, 0, 386),
    (5, 1, 380),
    (7, 0, 136),
)
TERMS = tuple(
    (iv.mpf(a) + b * ROOT3, coefficient)
    for a, b, coefficient in DISTANCE_COEFFICIENTS
)

# (first node, one past last node, interval upper bound)
NODE_BLOCKS = (
    (0, 5000, -iv.mpf(247) / 250),
    (5000, 10000, iv.mpf(-12)),
    (10000, 15000, -iv.mpf(99) / 100),
    (15000, 20000, iv.mpf(-18)),
    (20000, 25000, iv.mpf(-23)),
    (25000, 30000, iv.mpf(-22)),
    (30000, 35000, iv.mpf(-23)),
    (35000, 40001, iv.mpf(-26)),
)


def sinc_at_node(distance_squared, node):
    if node == 0:
        return iv.mpf(1)
    t = iv.mpf(node) / NODE_DENOMINATOR
    x = iv.sqrt(distance_squared) * t
    return iv.sin(x) / x


def verify_analytic_bounds():
    sqrt3_lower = Fraction(5, 3)
    sqrt3_upper = Fraction(26, 15)
    assert sqrt3_lower**2 < 3 < sqrt3_upper**2

    def upper_squared_distance(a, b):
        endpoint = sqrt3_upper if b >= 0 else sqrt3_lower
        return Fraction(a) + b * endpoint

    assert all(
        upper_squared_distance(a, b) < 9
        for a, b, _ in DISTANCE_COEFFICIENTS
    )
    coefficient_sum = sum(
        abs(coefficient) for _, _, coefficient in DISTANCE_COEFFICIENTS
    )
    assert coefficient_sum == 620268
    assert Fraction(3 * coefficient_sum, 2000) < 931

    assert 2 - sqrt3_upper > Fraction(1, 4)
    tail_upper = -Fraction(36608, 1000) + Fraction(2 * coefficient_sum, 40000)
    assert tail_upper == -Fraction(27973, 5000) < -Fraction(209, 400)


def main():
    assert NODE_BLOCKS[0][0] == 0
    assert NODE_BLOCKS[-1][1] == NODE_COUNT
    assert all(
        left[1] == right[0]
        for left, right in zip(NODE_BLOCKS, NODE_BLOCKS[1:])
    )
    for start, stop, bound in NODE_BLOCKS:
        for node in range(start, stop):
            value = iv.mpf(-36608)
            for distance_squared, coefficient in TERMS:
                value += coefficient * sinc_at_node(distance_squared, node)
            upper = (value / 1000).b
            assert upper < bound, (node, upper, bound)

    verify_analytic_bounds()


if __name__ == "__main__":
    main()
