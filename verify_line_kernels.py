#!/usr/bin/env python3
"""Verify the sampled higher-dimensional line kernel bounds.

The intervals use fixed-point integers.  Every returned endpoint is divided
by ``UNIT``; Python integers preserve the directed rounding without floating
point arithmetic.
"""

from dataclasses import dataclass
from decimal import Decimal, localcontext


DIGITS = 140
UNIT = 10**DIGITS
TOLERANCE = 10 ** (DIGITS - 35)


@dataclass(frozen=True)
class Check:
    name: str
    dimension: int
    grid_denominator: int
    last_node: int
    bound_numerator: int
    bound_denominator: int
    shifts: tuple[int, ...]
    coefficients: tuple[int, ...]


CHECKS = (
    Check("(ell_3,ell_4)", 3, 5000, 30000, 69, 100,
          (1, 2), (-3, -2)),
    Check("(ell_4,ell_4)", 4, 2500, 62500, 887, 1000,
          (1, 2, 3, 4), (-6, -4, -2, 1)),
    Check("(ell_3,ell_6)", 7, 10000, 100000, 399, 1000,
          (1, 2, 4, 5, 6), (-9, -5, -2, -1, 1)),
)


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def omega_interval(check: Check, numerator: int) -> tuple[int, int]:
    """Return a fixed-point interval for the relevant Omega function."""
    term_lower = term_upper = UNIT
    sum_lower = sum_upper = UNIT
    k = 0

    while True:
        n = k + 1
        denominator = (
            2 * check.grid_denominator**2 * n
            * (check.dimension + 2 * k)
        )
        next_lower = term_lower * numerator // denominator
        next_upper = ceil_div(term_upper * numerator, denominator)

        if n & 1:
            sum_lower -= next_upper
            sum_upper -= next_lower
        else:
            sum_lower += next_lower
            sum_upper += next_upper

        term_lower, term_upper = next_lower, next_upper
        k = n

        omitted_index = k + 1
        omitted_denominator = (
            2 * check.grid_denominator**2 * omitted_index
            * (check.dimension + 2 * k)
        )
        omitted_upper = ceil_div(
            term_upper * numerator, omitted_denominator
        )
        following_index = omitted_index + 1
        monotonicity_denominator = (
            2 * check.grid_denominator**2 * following_index
            * (check.dimension + 2 * omitted_index)
        )

        if (
            numerator <= monotonicity_denominator
            and omitted_upper < TOLERANCE
        ):
            if omitted_index & 1:
                sum_lower -= omitted_upper
            else:
                sum_upper += omitted_upper
            return sum_lower, sum_upper


def decimal_value(value: int) -> str:
    with localcontext() as context:
        context.prec = 80
        return format(Decimal(value) / Decimal(UNIT), ".60f")


def verify(check: Check) -> None:
    assert len(check.shifts) == len(check.coefficients)
    extremal_endpoint = -UNIT
    extremal_node = 0

    for node in range(check.last_node + 1):
        endpoint = 0
        node_squared = node * node
        for shift, coefficient in zip(check.shifts, check.coefficients):
            numerator = shift * shift * node_squared
            lower, upper = omega_interval(check, numerator)
            endpoint += coefficient * (upper if coefficient >= 0 else lower)

        if endpoint > extremal_endpoint:
            extremal_endpoint = endpoint
            extremal_node = node

        assert (
            endpoint * check.bound_denominator
            < UNIT * check.bound_numerator
        ), (check.name, node, endpoint)

    print(check.name)
    print(f"nodes verified: {check.last_node + 1}")
    print(f"extremal node: {extremal_node}")
    print(f"extremal upper endpoint: {decimal_value(extremal_endpoint)}")


def main() -> None:
    for check in CHECKS:
        verify(check)


if __name__ == "__main__":
    main()
