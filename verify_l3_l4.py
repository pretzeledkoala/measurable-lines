#!/usr/bin/env python3
"""Verify the finite and analytic bounds for ``(ell_3, ell_4)``."""

from fractions import Fraction

from _line_local import check


PHI = {
    (0, 0, 0): 0,
    (0, 0, 1): -1,
    (0, 1, 0): -2,
    (0, 1, 1): -2,
    (1, 0, 0): -3,
    (1, 0, 1): -3,
    (1, 1, 0): -3,
}


def allowed(word):
    return (
        (1, 1, 1) not in zip(word, word[1:], word[2:])
        and word != (0, 0, 0, 0)
    )


def local_weight(word):
    return (
        4 * word[0]
        - 1
        - 3 * word[0] * word[1]
        - 2 * word[0] * word[2]
    )


def verify_rational_bounds():
    compact = -Fraction(69, 100) - Fraction(7, 2) / 10000
    assert compact == -Fraction(13807, 20000) > -Fraction(7, 10)
    assert -Fraction(4, 6) > -Fraction(7, 10)

    square = -Fraction(57, 10)
    center = Fraction(47, 114)
    constant = -Fraction(71, 2280)
    assert -2 * square * center == Fraction(47, 10)
    assert square * center * center + constant == -1


def main():
    check((0, 1), 3, allowed, local_weight, PHI, (7, 12))
    verify_rational_bounds()


if __name__ == "__main__":
    main()
