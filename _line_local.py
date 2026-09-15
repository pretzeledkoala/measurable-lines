"""Shared exact checks for the higher-dimensional line inequalities."""

from itertools import product


def check(alphabet, length, allowed, weight, phi, expected):
    states = [
        word for word in product(alphabet, repeat=length) if allowed(word)
    ]
    state_set = set(states)
    transitions = [
        (state, state[1:] + (symbol,))
        for state in states
        for symbol in alphabet
        if allowed(state + (symbol,))
    ]
    assert state_set == set(phi)
    assert all(target in state_set for _, target in transitions)

    slacks = []
    for source, target in transitions:
        word = source + (target[-1],)
        slack = weight(word) + phi[source] - phi[target]
        assert slack >= 0, (word, slack)
        slacks.append(slack)

    assert (len(states), len(slacks)) == expected
    assert min(slacks) == 0
