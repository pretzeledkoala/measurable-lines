#!/usr/bin/env python3
"""Verify the finite square--tetrahedron inequality in Q(sqrt(3))."""

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


@dataclass(frozen=True)
class Q3:
    """The number a + b*sqrt(3), with a,b rational."""

    a: Fraction
    b: Fraction = Fraction(0)

    def __add__(self, other):
        return Q3(self.a + other.a, self.b + other.b)

    def __sub__(self, other):
        return Q3(self.a - other.a, self.b - other.b)

    def __neg__(self):
        return Q3(-self.a, -self.b)

    def __mul__(self, other):
        return Q3(
            self.a * other.a + 3 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )


ZERO = Q3(Fraction(0))
ONE = Q3(Fraction(1))
TWO = Q3(Fraction(2))
HALF = Q3(Fraction(1, 2))
ROOT3_OVER_2 = Q3(Fraction(0), Fraction(1, 2))
Point = tuple[Q3, Q3, Q3]


def point(x=ZERO, y=ZERO, z=ZERO):
    return (x, y, z)


def add(p: Point, q: Point) -> Point:
    return tuple(x + y for x, y in zip(p, q))


def squared_distance(p: Point, q: Point) -> Q3:
    out = ZERO
    for x, y in zip(p, q):
        d = x - y
        out = out + d * d
    return out


def make_points() -> list[Point]:
    e = point(ONE, ZERO, ZERO)
    v0 = point(ZERO, ONE, ZERO)
    w0 = point(ZERO, ZERO, ONE)
    v1 = point(ZERO, HALF, ROOT3_OVER_2)
    w1 = point(ZERO, -ROOT3_OVER_2, HALF)
    v2 = point(ZERO, -HALF, ROOT3_OVER_2)
    w2 = point(ZERO, -ROOT3_OVER_2, -HALF)

    return [
        point(), w0, v0, add(v0, w0),
        e, add(e, w0), add(e, v0), add(e, add(v0, w0)),
        w1, v1, add(v1, w1), add(e, w1), add(e, v1),
        add(e, add(v1, w1)), w2, v2, add(v2, w2), add(e, w2),
        add(e, v2), add(e, add(v2, w2)),
    ]


def permutation_from_cycles(
    n: int, cycles: list[tuple[int, ...]]
) -> tuple[int, ...]:
    p = list(range(n))
    for cycle in cycles:
        for x, y in zip(cycle, cycle[1:] + cycle[:1]):
            p[x] = y
    return tuple(p)


def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(p[q[i]] for i in range(len(p)))


def make_group() -> list[tuple[int, ...]]:
    sigma = permutation_from_cycles(20, [
        (1, 15), (2, 14), (3, 16), (5, 18), (6, 17), (7, 19),
        (8, 9), (11, 12),
    ])
    tau = permutation_from_cycles(20, [
        (0, 4), (1, 5), (2, 6), (3, 7), (8, 11), (9, 12),
        (10, 13), (14, 17), (15, 18), (16, 19),
    ])

    identity = tuple(range(20))
    group = {identity}
    changed = True
    while changed:
        changed = False
        for g in tuple(group):
            for h in (sigma, tau):
                for gh in (compose(g, h), compose(h, g)):
                    if gh not in group:
                        group.add(gh)
                        changed = True
    return sorted(group)


# representative, orbit size, squared distance, coefficient
PAIR_ORBIT_DATA = [
    ((8, 10), 4, Q3(Fraction(1)), 0),
    ((0, 8), 4, Q3(Fraction(1)), -17344),
    ((0, 10), 2, Q3(Fraction(2)), -152),
    ((0, 1), 4, Q3(Fraction(1)), -9839),
    ((4, 8), 4, Q3(Fraction(2)), -21087),
    ((0, 2), 4, Q3(Fraction(1)), -10680),
    ((4, 10), 2, Q3(Fraction(3)), 587),
    ((0, 3), 4, Q3(Fraction(2)), -2319),
    ((3, 4), 4, Q3(Fraction(3)), 2040),
    ((1, 4), 4, Q3(Fraction(2)), -553),
    ((1, 3), 4, Q3(Fraction(1)), 3317),
    ((2, 3), 4, Q3(Fraction(1)), -997),
    ((7, 9), 4, Q3(Fraction(3), Fraction(-1)), 5100),
    ((0, 4), 1, Q3(Fraction(1)), -39099),
    ((5, 14), 4, Q3(Fraction(4)), -1425),
    ((2, 4), 4, Q3(Fraction(2)), -1277),
    ((2, 14), 2, Q3(Fraction(2), Fraction(1)), -12),
    ((1, 14), 4, Q3(Fraction(3)), 1304),
    ((3, 14), 4, Q3(Fraction(4), Fraction(1)), 51),
    ((6, 14), 2, Q3(Fraction(3), Fraction(1)), -22),
    ((1, 9), 4, Q3(Fraction(2), Fraction(-1)), -6903),
    ((7, 14), 4, Q3(Fraction(5), Fraction(1)), 95),
    ((1, 2), 4, Q3(Fraction(2)), -1171),
    ((5, 9), 4, Q3(Fraction(3), Fraction(-1)), 5627),
    ((1, 5), 2, Q3(Fraction(1)), -6775),
    ((5, 10), 4, Q3(Fraction(3), Fraction(-1)), 485),
    ((2, 5), 4, Q3(Fraction(3)), 1144),
    ((2, 9), 4, Q3(Fraction(1)), -9068),
    ((1, 8), 4, Q3(Fraction(1)), -2603),
    ((3, 5), 4, Q3(Fraction(2)), -2761),
    ((3, 15), 4, Q3(Fraction(4), Fraction(-1)), -3658),
    ((1, 10), 4, Q3(Fraction(2), Fraction(-1)), -337),
    ((3, 16), 2, Q3(Fraction(6)), 193),
    ((7, 15), 4, Q3(Fraction(5), Fraction(-1)), 2749),
    ((7, 16), 2, Q3(Fraction(7)), 68),
    ((2, 10), 4, Q3(Fraction(2), Fraction(1)), -294),
    ((3, 9), 4, Q3(Fraction(2), Fraction(-1)), -6010),
    ((3, 7), 2, Q3(Fraction(1)), -2480),
    ((9, 11), 2, Q3(Fraction(3)), 460),
    ((5, 8), 4, Q3(Fraction(2)), -5565),
    ((1, 15), 2, Q3(Fraction(2), Fraction(-1)), -7650),
    ((6, 8), 4, Q3(Fraction(3), Fraction(1)), -21),
    ((7, 8), 4, Q3(Fraction(3), Fraction(1)), -145),
    ((10, 11), 4, Q3(Fraction(2)), 537),
    ((2, 8), 4, Q3(Fraction(2), Fraction(1)), -93),
    ((3, 10), 4, Q3(Fraction(2)), -59),
    ((3, 6), 4, Q3(Fraction(2)), 1433),
    ((6, 9), 4, Q3(Fraction(2)), -2065),
    ((3, 8), 4, Q3(Fraction(2), Fraction(1)), 550),
    ((5, 15), 2, Q3(Fraction(3), Fraction(-1)), 5870),
    ((6, 10), 4, Q3(Fraction(3), Fraction(1)), 76),
    ((8, 9), 2, Q3(Fraction(2)), -1510),
    ((7, 10), 4, Q3(Fraction(3)), -428),
    ((10, 13), 1, Q3(Fraction(1)), -267),
    ((2, 6), 2, Q3(Fraction(1)), -54),
    ((8, 11), 2, Q3(Fraction(1)), -26934),
]
PAPER_LOCAL_NUMERATOR_BOUND = 0

EXPECTED_DISTANCE_COEFFICIENTS = {
    Q3(Fraction(2), Fraction(-1)): -68300,
    Q3(Fraction(1)): -300708,
    Q3(Fraction(3), Fraction(-1)): 56588,
    Q3(Fraction(2)): -142872,
    Q3(Fraction(4), Fraction(-1)): -14632,
    Q3(Fraction(3)): 18334,
    Q3(Fraction(5), Fraction(-1)): 10996,
    Q3(Fraction(2), Fraction(1)): 628,
    Q3(Fraction(4)): -5700,
    Q3(Fraction(3), Fraction(1)): -404,
    Q3(Fraction(4), Fraction(1)): 204,
    Q3(Fraction(6)): 386,
    Q3(Fraction(5), Fraction(1)): 380,
    Q3(Fraction(7)): 136,
}


def pair_orbit(rep: tuple[int, int], group: list[tuple[int, ...]]):
    return frozenset(tuple(sorted((g[rep[0]], g[rep[1]]))) for g in group)


def signed_pair_sum(mask: int, orbit) -> int:
    total = 0
    for i, j in orbit:
        total += 1 if ((mask >> i) & 1) == ((mask >> j) & 1) else -1
    return total


def classify_candidate(points: list[Point], subset: tuple[int, int, int, int]):
    """Classify a four-set with distances in {1,sqrt(2)}.

    Returns ``"square"`` when the four unit edges form a 4-cycle and the
    two sqrt(2) edges are disjoint diagonals.  Returns ``"tetrahedron"``
    when the two sqrt(2) edges are adjacent.  Otherwise returns
    ``None``.
    """
    pairs = list(combinations(subset, 2))
    unit_edges = []
    diagonal_edges = []
    for i, j in pairs:
        distance = squared_distance(points[i], points[j])
        if distance == ONE:
            unit_edges.append((i, j))
        elif distance == TWO:
            diagonal_edges.append((i, j))
        else:
            return None

    if len(unit_edges) != 4 or len(diagonal_edges) != 2:
        return None

    if set(diagonal_edges[0]).isdisjoint(diagonal_edges[1]):
        degrees = {vertex: 0 for vertex in subset}
        for i, j in unit_edges:
            degrees[i] += 1
            degrees[j] += 1
        if all(degree == 2 for degree in degrees.values()):
            return "square"
        return None

    return "tetrahedron"


def main() -> None:
    points = make_points()
    assert len(points) == 20
    assert len(set(points)) == 20

    squares = []
    tetrahedra = []
    for subset in combinations(range(20), 4):
        kind = classify_candidate(points, subset)
        if kind == "square":
            squares.append(subset)
        elif kind == "tetrahedron":
            tetrahedra.append(subset)

    assert len(squares) == 22
    assert len(tetrahedra) == 24
    forbidden_configurations = squares + tetrahedra
    assert len(forbidden_configurations) == 46

    group = make_group()
    assert len(group) == 4
    for g in group:
        for i in range(20):
            for j in range(20):
                assert (
                    squared_distance(points[i], points[j])
                    == squared_distance(points[g[i]], points[g[j]])
                )

    orbits = []
    for representative, expected_size, expected_distance, _ in PAIR_ORBIT_DATA:
        orbit = pair_orbit(representative, group)
        assert len(orbit) == expected_size
        actual_distances = {
            squared_distance(points[i], points[j]) for i, j in orbit
        }
        assert actual_distances == {expected_distance}
        orbits.append(orbit)

    assert len(set(orbits)) == 56
    assert sum(map(len, orbits)) == 190

    grouped = {}
    for orbit, (_, _, distance, coefficient) in zip(orbits, PAIR_ORBIT_DATA):
        grouped[distance] = grouped.get(distance, 0) + len(orbit) * coefficient
    assert grouped == EXPECTED_DISTANCE_COEFFICIENTS
    assert sum(abs(value) for value in grouped.values()) == 620268

    forbidden_masks = [
        sum(1 << i for i in configuration)
        for configuration in forbidden_configurations
    ]
    allowed_count = 0
    minimum = None
    minimizers = []

    for mask in range(1 << 20):
        if any((mask & configuration) in (0, configuration)
               for configuration in forbidden_masks):
            continue

        allowed_count += 1
        numerator = -36608
        for orbit, (_, _, _, coefficient) in zip(orbits, PAIR_ORBIT_DATA):
            numerator += coefficient * signed_pair_sum(mask, orbit)
        assert numerator >= PAPER_LOCAL_NUMERATOR_BOUND, (mask, numerator)

        if minimum is None or numerator < minimum:
            minimum = numerator
            minimizers = [mask]
        elif numerator == minimum:
            minimizers.append(mask)

    assert allowed_count == 47946
    assert minimum == PAPER_LOCAL_NUMERATOR_BOUND
    assert len(minimizers) == 4


if __name__ == "__main__":
    main()
