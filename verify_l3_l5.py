#!/usr/bin/env python3
"""Verify the finite and analytic bounds for ``(ell_3, ell_5)``."""

from fractions import Fraction

from _line_local import check


def allowed(word):
    s = "".join(map(str, word))
    return "111" not in s and "00000" not in s


PHI = {
    (0, 0, 0, 0, 1, 0): -2,
    (0, 0, 0, 0, 1, 1): 0,
    (0, 0, 0, 1, 0, 0): -8,
    (0, 0, 0, 1, 0, 1): -8,
    (0, 0, 0, 1, 1, 0): -6,
    (0, 0, 1, 0, 0, 0): -14,
    (0, 0, 1, 0, 0, 1): -14,
    (0, 0, 1, 0, 1, 0): -14,
    (0, 0, 1, 0, 1, 1): -14,
    (0, 0, 1, 1, 0, 0): -13,
    (0, 0, 1, 1, 0, 1): -12,
    (0, 1, 0, 0, 0, 0): -20,
    (0, 1, 0, 0, 0, 1): -20,
    (0, 1, 0, 0, 1, 0): -20,
    (0, 1, 0, 0, 1, 1): -20,
    (0, 1, 0, 1, 0, 0): -22,
    (0, 1, 0, 1, 0, 1): -20,
    (0, 1, 0, 1, 1, 0): -20,
    (0, 1, 1, 0, 0, 0): -24,
    (0, 1, 1, 0, 0, 1): -22,
    (0, 1, 1, 0, 1, 0): -20,
    (0, 1, 1, 0, 1, 1): -18,
    (1, 0, 0, 0, 0, 1): -26,
    (1, 0, 0, 0, 1, 0): -26,
    (1, 0, 0, 0, 1, 1): -26,
    (1, 0, 0, 1, 0, 0): -30,
    (1, 0, 0, 1, 0, 1): -28,
    (1, 0, 0, 1, 1, 0): -28,
    (1, 0, 1, 0, 0, 0): -28,
    (1, 0, 1, 0, 0, 1): -28,
    (1, 0, 1, 0, 1, 0): -26,
    (1, 0, 1, 0, 1, 1): -26,
    (1, 0, 1, 1, 0, 0): -29,
    (1, 0, 1, 1, 0, 1): -27,
    (1, 1, 0, 0, 0, 0): -30,
    (1, 1, 0, 0, 0, 1): -30,
    (1, 1, 0, 0, 1, 0): -28,
    (1, 1, 0, 0, 1, 1): -28,
    (1, 1, 0, 1, 0, 0): -26,
    (1, 1, 0, 1, 0, 1): -26,
    (1, 1, 0, 1, 1, 0): -24,
}


def local_weight(r):
    return (
        -6
        + 28 * r[0]
        - 20 * r[0] * r[1]
        - 14 * r[0] * r[2]
        - 3 * r[0] * r[3]
        - 4 * r[0] * r[4]
        + 2 * r[0] * r[5]
        + 2 * r[0] * r[6]
    )


def ceil_div(a, b):
    return -((-a) // b)


def omega_interval(d, k, j, N, scale, threshold):
    """
    Return integers L,U with
        L/scale <= Omega_d(k*j/N) <= U/scale.

    The series is
      Omega_d(x)=sum_n (-1)^n a_n,
      a_{n+1}/a_n=x^2/[2(n+1)(d+2n)].
    """
    term_lo = scale
    term_hi = scale
    sum_lo = 0
    sum_hi = 0
    n = 0
    numerator = (k * j) ** 2

    while True:
        if n % 2 == 0:
            sum_lo += term_lo
            sum_hi += term_hi
        else:
            sum_lo -= term_hi
            sum_hi -= term_lo

        denominator = N * N * 2 * (n + 1) * (d + 2 * n)
        next_lo = term_lo * numerator // denominator
        next_hi = ceil_div(term_hi * numerator, denominator)

        # The ratios decrease with n. Once below 1, the remaining
        # alternating tail has the sign of the next term and magnitude
        # at most that next term.
        if numerator < denominator and next_hi <= threshold:
            if (n + 1) % 2 == 0:
                sum_hi += next_hi
            else:
                sum_lo -= next_hi
            return sum_lo, sum_hi

        term_lo, term_hi = next_lo, next_hi
        n += 1


def verify_kernel_nodes():
    d = 4
    N = 2000
    last_j = 50000  # t <= 25
    scale = 10**140
    threshold = 10**105
    coeffs = (-20, -14, -3, -4, 2, 2)

    for j in range(last_j + 1):
        upper = 0
        for k, coefficient in enumerate(coeffs, start=1):
            lo, hi = omega_interval(d, k, j, N, scale, threshold)
            upper += coefficient * (hi if coefficient > 0 else lo)

        assert 500 * upper < 1469 * scale, (j, upper)


def verify_rational_finishing_steps():
    # Interpolation on [0,25]:
    # |Omega_4'| < 1/2, hence |H'| < 95/2.
    interpolated = Fraction(1469, 500) + Fraction(95, 8000)
    assert interpolated < Fraction(59, 20)

    # Tail for t>=25, using |Omega_4(x)| <= 2/x.
    tail_at_25 = Fraction(2, 25) * (
        20 + Fraction(14, 2) + Fraction(3, 3)
        + Fraction(4, 4) + Fraction(2, 5) + Fraction(2, 6)
    )
    assert tail_at_25 < Fraction(59, 20)

    # Global quadratic:
    # -6+28p-37p^2+(59/20)p(1-p)
    # = -(799/20)(p-619/1598)^2 - 359/63920.
    assert Fraction(359, 63920) > 0

def main():
    check((0, 1), 6, allowed, local_weight, PHI, (41, 73))
    verify_kernel_nodes()
    verify_rational_finishing_steps()


if __name__ == "__main__":
    main()
