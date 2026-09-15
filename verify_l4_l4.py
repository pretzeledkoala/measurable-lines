#!/usr/bin/env python3
"""Verify the finite and analytic bounds for ``(ell_4, ell_4)``."""

from fractions import Fraction

from _line_local import check


PHI = {
    (-1, -1, -1, 1): 0,
    (-1, -1, 1, -1): -8,
    (-1, -1, 1, 1): -10,
    (-1, 1, -1, -1): -12,
    (-1, 1, -1, 1): -14,
    (-1, 1, 1, -1): -14,
    (-1, 1, 1, 1): -12,
    (1, -1, -1, -1): -12,
    (1, -1, -1, 1): -14,
    (1, -1, 1, -1): -14,
    (1, -1, 1, 1): -12,
    (1, 1, -1, -1): -10,
    (1, 1, -1, 1): -8,
    (1, 1, 1, -1): 0,
}


def allowed(word):
    last = word[-4:]
    return not (all(value == 1 for value in last) or
                all(value == -1 for value in last))


def local_weight(word):
    return (
        -1
        - 6 * word[0] * word[1]
        - 4 * word[0] * word[2]
        - 2 * word[0] * word[3]
        + word[0] * word[4]
    )


def verify_rational_bounds():
    compact = Fraction(887, 1000) + Fraction(12, 5000)
    assert compact == Fraction(4447, 5000) < Fraction(9, 10)
    tail = Fraction(107, 6 * 25)
    assert tail == Fraction(107, 150) < Fraction(9, 10)


def main():
    check((-1, 1), 4, allowed, local_weight, PHI, (14, 26))
    verify_rational_bounds()


if __name__ == "__main__":
    main()
