#!/usr/bin/env python3
"""Verify the finite and analytic bounds for ``(ell_3, ell_6)``."""

from fractions import Fraction

from _line_local import check


PHI = {
    (0, 0, 0, 0, 0, 1): 0,
    (0, 0, 0, 0, 1, 0): -2,
    (0, 0, 0, 0, 1, 1): -2,
    (0, 0, 0, 1, 0, 0): -4,
    (0, 0, 0, 1, 0, 1): -4,
    (0, 0, 0, 1, 1, 0): -4,
    (0, 0, 1, 0, 0, 0): -6,
    (0, 0, 1, 0, 0, 1): -6,
    (0, 0, 1, 0, 1, 0): -6,
    (0, 0, 1, 0, 1, 1): -6,
    (0, 0, 1, 1, 0, 0): -6,
    (0, 0, 1, 1, 0, 1): -6,
    (0, 1, 0, 0, 0, 0): -8,
    (0, 1, 0, 0, 0, 1): -8,
    (0, 1, 0, 0, 1, 0): -8,
    (0, 1, 0, 0, 1, 1): -8,
    (0, 1, 0, 1, 0, 0): -9,
    (0, 1, 0, 1, 0, 1): -8,
    (0, 1, 0, 1, 1, 0): -9,
    (0, 1, 1, 0, 0, 0): -8,
    (0, 1, 1, 0, 0, 1): -8,
    (0, 1, 1, 0, 1, 0): -8,
    (0, 1, 1, 0, 1, 1): -8,
    (1, 0, 0, 0, 0, 0): -10,
    (1, 0, 0, 0, 0, 1): -10,
    (1, 0, 0, 0, 1, 0): -11,
    (1, 0, 0, 0, 1, 1): -10,
    (1, 0, 0, 1, 0, 0): -12,
    (1, 0, 0, 1, 0, 1): -11,
    (1, 0, 0, 1, 1, 0): -13,
    (1, 0, 1, 0, 0, 0): -11,
    (1, 0, 1, 0, 0, 1): -11,
    (1, 0, 1, 0, 1, 0): -11,
    (1, 0, 1, 0, 1, 1): -10,
    (1, 0, 1, 1, 0, 0): -12,
    (1, 0, 1, 1, 0, 1): -11,
    (1, 1, 0, 0, 0, 0): -10,
    (1, 1, 0, 0, 0, 1): -10,
    (1, 1, 0, 0, 1, 0): -10,
    (1, 1, 0, 0, 1, 1): -10,
    (1, 1, 0, 1, 0, 0): -10,
    (1, 1, 0, 1, 0, 1): -10,
    (1, 1, 0, 1, 1, 0): -10,
}


def allowed(word):
    return (
        (1, 1, 1) not in zip(word, word[1:], word[2:])
        and word[-6:] != (0, 0, 0, 0, 0, 0)
    )


def local_weight(word):
    return (
        -2
        + 11 * word[0]
        - 9 * word[0] * word[1]
        - 5 * word[0] * word[2]
        - 2 * word[0] * word[4]
        - word[0] * word[5]
        + word[0] * word[6]
    )


def verify_rational_bounds():
    compact = Fraction(399, 1000) + Fraction(95, 8) / 20000
    assert compact == Fraction(12787, 32000) < Fraction(2, 5)
    tail = Fraction(30, 10**3) * (
        9
        + Fraction(5, 2**3)
        + Fraction(2, 4**3)
        + Fraction(1, 5**3)
        + Fraction(1, 6**3)
    )
    assert tail == Fraction(1044239, 3600000) < Fraction(2, 5)

    square = -Fraction(82, 5)
    center = Fraction(57, 164)
    constant = -Fraction(31, 1640)
    assert -2 * square * center == Fraction(57, 5)
    assert square * center * center + constant == -2


def main():
    check((0, 1), 6, allowed, local_weight, PHI, (43, 78))
    verify_rational_bounds()


if __name__ == "__main__":
    main()
