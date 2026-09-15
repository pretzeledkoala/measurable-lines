#!/usr/bin/env python3
"""Verify the bounded-interval kernel estimate for ``(ell_2, ell_7)``."""

from fractions import Fraction


NK = 14
DENOMINATOR = 1_000_000
TMAX = 100
SCALE = 10**420
REM_LIMIT = 10**340

KERNEL_Q = (1, 3, 4, 7, 9, 13, 16, 21, 25, 31, 36, 43, 49, 57)
KERNEL_COEFFICIENTS = (
    -489, 1035, -669, 600, -700, 548, -355,
    441, -105, 339, -114, -139, 16, -90,
)
KERNEL_CONSTANT = -319

# |K'(t)| <= DERIVATIVE_NUMERATOR / DERIVATIVE_DENOMINATOR.
DERIVATIVE_NUMERATOR = 35046648006
DERIVATIVE_DENOMINATOR = 3000000

SQRT_UPPER_NUMERATOR = (
    1000000, 1732051, 2000000, 2645752, 3000000, 3605552,
    4000000, 4582576, 5000000, 5567765, 6000000, 6557439,
    7000000, 7549835,
)
FOURTH_LOWER_NUMERATOR = (
    1, 13, 7, 8, 17, 189, 2, 21, 11, 47, 12, 5, 13, 27,
)
FOURTH_LOWER_DENOMINATOR = (
    1, 10, 5, 5, 10, 100, 1, 10, 5, 20, 5, 2, 5, 10,
)


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def verify_rational_bounds() -> None:
    derivative_sum = 0
    tail = Fraction(0)
    coefficient_sum = KERNEL_CONSTANT

    for q, coefficient, sqrt_upper, fourth_lower, denominator in zip(
        KERNEL_Q,
        KERNEL_COEFFICIENTS,
        SQRT_UPPER_NUMERATOR,
        FOURTH_LOWER_NUMERATOR,
        FOURTH_LOWER_DENOMINATOR,
    ):
        assert sqrt_upper**2 >= q * 10**12
        derivative_sum += abs(coefficient) * sqrt_upper

        assert fourth_lower**4 <= q * denominator**4
        tail += Fraction(abs(coefficient) * denominator, fourth_lower)
        coefficient_sum += coefficient

    assert coefficient_sum == -1
    assert 2 * derivative_sum == DERIVATIVE_NUMERATOR
    assert tail < 3600
    print(f"tail coefficient sum upper={tail} < 3600")


def j0_interval(q: int, node: int) -> tuple[int, int]:
    """Return [lower, upper] scaled by SCALE for J_0(sqrt(q)t)."""
    if node == 0:
        return SCALE, SCALE

    numerator = q * node * node
    denominator_base = 4 * DENOMINATOR**2
    term_lower = term_upper = SCALE
    sum_lower = sum_upper = SCALE
    k = 0

    while True:
        next_index = k + 1
        denominator = denominator_base * next_index**2
        next_lower = term_lower * numerator // denominator
        next_upper = ceil_div(term_upper * numerator, denominator)

        omitted_index = next_index + 1
        decreasing_denominator = denominator_base * omitted_index**2
        if (
            numerator <= decreasing_denominator
            and next_upper < REM_LIMIT
        ):
            if next_index & 1:
                sum_lower -= next_upper
            else:
                sum_upper += next_upper
            return sum_lower, sum_upper

        k = next_index
        term_lower, term_upper = next_lower, next_upper
        if k & 1:
            sum_lower -= term_upper
            sum_upper -= term_lower
        else:
            sum_lower += term_lower
            sum_upper += term_upper
        assert k < 2000


def kernel_upper(node: int) -> int:
    upper = SCALE * KERNEL_CONSTANT
    for q, coefficient in zip(KERNEL_Q, KERNEL_COEFFICIENTS):
        lower, upper_term = j0_interval(q, node)
        upper += coefficient * (upper_term if coefficient > 0 else lower)
    return upper


def verify_kernel() -> None:
    end = TMAX * DENOMINATOR
    node = 0
    nodes = 0
    minimum_step = None
    maximum_step = 0
    step_denominator = 2 * SCALE * DERIVATIVE_NUMERATOR

    while node < end:
        upper = kernel_upper(node)
        assert upper < 0, (node, upper)

        step = (
            -upper * DENOMINATOR * DERIVATIVE_DENOMINATOR
            // step_denominator
        )
        assert step > 0, (node, upper)
        step = min(step, DENOMINATOR // 10, end - node)

        minimum_step = step if minimum_step is None else min(minimum_step, step)
        maximum_step = max(maximum_step, step)
        node += step
        nodes += 1

        if nodes % 1000 == 0:
            print(f"kernel nodes={nodes} t={node / DENOMINATOR:.6f}")

    assert kernel_upper(end) < 0
    print(
        f"kernel verified on [0,{TMAX}], nodes={nodes}, "
        f"step range=[{minimum_step},{maximum_step}]/{DENOMINATOR}"
    )


def main() -> None:
    assert len(KERNEL_Q) == len(KERNEL_COEFFICIENTS) == NK
    verify_rational_bounds()
    verify_kernel()
    print("all checks passed")


if __name__ == "__main__":
    main()
