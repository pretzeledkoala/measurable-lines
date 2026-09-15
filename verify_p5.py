#!/usr/bin/env python3
"""Verify the finite and kernel data for P_5.

The local inequalities use exact arithmetic.  The compact kernel checks use
rational interval bounds for the alternating Bessel series; the analytic
tail estimates are checked separately below.
"""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


F = Fraction


@dataclass(frozen=True)
class Quadratic:
    """An element a+b*sqrt(d), with rational a and b."""

    a: Fraction
    b: Fraction
    d: int

    def _check(self, other):
        assert self.d == other.d

    def __add__(self, other):
        self._check(other)
        return Quadratic(self.a + other.a, self.b + other.b, self.d)

    def __mul__(self, other):
        self._check(other)
        return Quadratic(
            self.a * other.a + self.d * self.b * other.b,
            self.a * other.b + self.b * other.a,
            self.d,
        )

    def __rmul__(self, scalar):
        if isinstance(scalar, Quadratic):
            return scalar * self
        scalar = F(scalar)
        return Quadratic(scalar * self.a, scalar * self.b, self.d)

def q5(a=0, b=0):
    return Quadratic(F(a), F(b), 5)


Q5_ZERO = q5()
Q5_ONE = q5(1)


def vector_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vector_sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def q5_dot_regular_vertices(i, j):
    # R_5^2=(5+sqrt(5))/10 and the three possible cosine values.
    radius_squared = q5(F(1, 2), F(1, 10))
    difference = (i - j) % 5
    difference = min(difference, 5 - difference)
    cosine = {
        0: q5(1),
        1: q5(F(-1, 4), F(1, 4)),
        2: q5(F(-1, 4), F(-1, 4)),
    }[difference]
    return radius_squared * cosine


def dot5(a, b):
    value = Q5_ZERO
    for i, coefficient_i in enumerate(a):
        for j, coefficient_j in enumerate(b):
            value += (
                coefficient_i
                * coefficient_j
                * q5_dot_regular_vertices(i, j)
            )
    return value


def unit_vector(index):
    return tuple(int(i == index) for i in range(5))


def p5_points():
    shift = vector_sub(unit_vector(1), unit_vector(0))
    translated = [vector_add(unit_vector(j), shift) for j in range(5)]
    assert translated[0] == unit_vector(1)
    return [unit_vector(j) for j in range(5)] + translated[1:]


def p5_squared_distance(a, b):
    difference = vector_sub(a, b)
    return dot5(difference, difference)


def allowed(mask, edges, blocks):
    if any(
        (mask >> i) & 1 and (mask >> j) & 1
        for i, j in edges
    ):
        return False
    return all(any((mask >> i) & 1 for i in block) for block in blocks)


def functional(mask, constant, linear, quadratic):
    value = constant
    value += sum(
        coefficient
        for i, coefficient in enumerate(linear)
        if (mask >> i) & 1
    )
    value += sum(
        coefficient
        for (i, j), coefficient in quadratic.items()
        if (mask >> i) & 1 and (mask >> j) & 1
    )
    return value


def bessel_interval(q_lower, q_upper, t, degree):
    """Enclose J_0(sqrt(q)*t) using rational interval arithmetic."""
    z_lower = q_lower * t * t / 4
    z_upper = q_upper * t * t / 4
    lower = upper = F(1)
    term_lower = term_upper = F(1)

    for m in range(1, degree + 1):
        denominator = m * m
        term_lower = term_lower * z_lower / denominator
        term_upper = term_upper * z_upper / denominator
        if m % 2:
            lower -= term_upper
            upper -= term_lower
        else:
            lower += term_lower
            upper += term_upper

    first_omitted = degree + 1
    omitted_lower = term_lower * z_lower / (first_omitted * first_omitted)
    omitted_upper = term_upper * z_upper / (first_omitted * first_omitted)
    ratio = z_upper / ((degree + 2) * (degree + 2))
    assert ratio < 1
    tail = omitted_upper / (1 - ratio)
    return lower - tail, upper + tail


def kernel_interval(constant, terms, t, degree):
    """Return a rational lower and upper bound for a Bessel kernel."""
    lower = upper = F(constant)
    for coefficient, q_lower, q_upper in terms:
        bessel_lower, bessel_upper = bessel_interval(
            q_lower, q_upper, t, degree
        )
        if coefficient >= 0:
            lower += coefficient * bessel_lower
            upper += coefficient * bessel_upper
        else:
            lower += coefficient * bessel_upper
            upper += coefficient * bessel_lower
    return lower, upper


def q5_interval(value, sqrt5_lower, sqrt5_upper):
    assert value.d == 5
    if value.b >= 0:
        return (
            value.a + value.b * sqrt5_lower,
            value.a + value.b * sqrt5_upper,
        )
    return (
        value.a + value.b * sqrt5_upper,
        value.a + value.b * sqrt5_lower,
    )


def verify_p5():
    p5_local_constant = -77
    points = p5_points()
    distances = {
        pair: p5_squared_distance(points[pair[0]], points[pair[1]])
        for pair in combinations(range(9), 2)
    }
    edges = tuple(pair for pair, distance in distances.items() if distance == Q5_ONE)
    expected_edges = (
        (0, 1), (0, 4), (1, 2), (1, 5), (1, 8), (2, 3),
        (2, 6), (3, 4), (3, 7), (4, 8), (5, 6), (6, 7), (7, 8),
    )
    assert edges == expected_edges
    blocks = ((0, 1, 2, 3, 4), (1, 5, 6, 7, 8))

    linear = (-23, 77, -100, -23, -23, 100, 100, 100, 77)
    rows = (
        (-100, 23, 23, -100, 0, 0, 0, 23),
        (-100, 23, 23, -100, -100, -100, -100),
        (-100, 23, 77, -100, 77, 100),
        (-100, 0, 0, -100, 23),
        (0, 0, 0, -100),
        (-100, -100, -100),
        (-48, -100),
        (100,),
    )
    quadratic = {}
    for i, row in enumerate(rows):
        for offset, coefficient in enumerate(row, start=i + 1):
            quadratic[(i, offset)] = coefficient

    valid = [
        mask for mask in range(1 << 9)
        if allowed(mask, edges, blocks)
    ]
    values_by_size = {}
    for mask in valid:
        size = mask.bit_count()
        value = functional(mask, p5_local_constant, linear, quadratic)
        count, minimum, maximum = values_by_size.get(size, (0, None, None))
        values_by_size[size] = (
            count + 1,
            value if minimum is None else min(minimum, value),
            value if maximum is None else max(maximum, value),
        )
    assert len(valid) == 42
    assert values_by_size == {
        1: (1, 0, 0),
        2: (17, 0, 0),
        3: (19, 0, 77),
        4: (5, 23, 100),
    }

    aggregates = Counter()
    for pair, coefficient in quadratic.items():
        if coefficient:
            aggregates[distances[pair]] += coefficient
    assert aggregates == {
        q5(F(3, 2), F(-1, 2)): 200,
        q5(1): -1048,
        q5(F(5, 2), F(-1, 2)): 100,
        q5(F(3, 2), F(1, 2)): -385,
    }
    assert sum(linear) == 285
    assert sum(quadratic.values()) == -1133
    assert sum(abs(c) for c in linear) == 623
    assert sum(abs(c) for c in quadratic.values()) == 2163
    assert 8 * 623 + 24 * 2163 == 56896

    sqrt5_lower = F(22360679774, 10**10)
    sqrt5_upper = F(22360679775, 10**10)
    assert sqrt5_lower * sqrt5_lower < 5 < sqrt5_upper * sqrt5_upper
    q_values = (
        q5(F(3, 2), F(-1, 2)), q5(1),
        q5(F(5, 2), F(-1, 2)), q5(F(3, 2), F(1, 2)),
    )
    q_intervals = tuple(
        q5_interval(q, sqrt5_lower, sqrt5_upper) for q in q_values
    )
    expected_q_intervals = (
        (F(3, 2) - sqrt5_upper / 2, F(3, 2) - sqrt5_lower / 2),
        (F(1), F(1)),
        (F(5, 2) - sqrt5_upper / 2, F(5, 2) - sqrt5_lower / 2),
        (F(3, 2) + sqrt5_lower / 2, F(3, 2) + sqrt5_upper / 2),
    )
    assert q_intervals == expected_q_intervals

    fourth_lower = (F(39, 50), F(1), F(27, 25), F(127, 100))
    assert all(
        lower**4 <= q_lower
        for (q_lower, _), lower in zip(q_intervals, fourth_lower)
    )
    assert all(
        lower**4 < q_lower
        for (q_lower, _), lower in zip(q_intervals, fourth_lower)
        if lower != 1
    )
    assert q_intervals[0][1] < 1
    assert q_intervals[2][1] < F(6, 5) ** 2
    assert q_intervals[3][1] < F(13, 8) ** 2

    terms = tuple(
        (coefficient, q_lower, q_upper)
        for coefficient, (q_lower, q_upper) in zip(
            (200, -1048, 100, -385), q_intervals
        )
    )
    derivative_sum = (
        200 + 1048 + 100 * F(6, 5) + 385 * F(13, 8)
    )
    assert F(2, 3) * derivative_sum == F(15949, 12)
    assert F(389) + F(15949, 4800) < 393
    tail_sum = (
        F(200) / F(39, 50) + 1048
        + F(100) / F(27, 25) + F(385) / F(127, 100)
    )
    assert tail_sum == F(75787696, 44577)
    assert F(25, 471) * tail_sum * tail_sum < 392**2

    maximum = None
    maximum_t = None
    for k in range(2400):
        t = F(2 * k + 1, 400)
        _, upper = kernel_interval(0, terms, t, 45)
        assert upper < 389, (t, upper)
        if maximum is None or upper > maximum:
            maximum, maximum_t = upper, t
    assert maximum < 389
    assert maximum_t is not None
    print("P5: 42 local states; 2400 nodes; every node K_5 < 389")
    print("P5: sqrt(5) in (22360679774/10^10, 22360679775/10^10)")


def main():
    verify_p5()


if __name__ == "__main__":
    main()
