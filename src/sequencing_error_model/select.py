"""Greedy criterion-based component selection for either head (plan §5.4, §9 phase 2).

Forward selection as in the fork's `model_selection.py`: screen the variants of the head's required first
component (`QualityMarkov(m)` for head Q, `QualityWindow(m)` for head E), then repeatedly add the candidate
group whose best variant lowers the criterion most, and stop when none does. Every criterion is scored on
a held-out table, which the caller builds from reads not used for training (split by read, not by base).

The default, `test-ll` (negative held-out log-likelihood), already guards against overfitting. `aic` and `bic`
add a penalty on top using raw parameter counts, which overcount (softmax gauge, masked cells, L2 shrinkage):
in the recovery tests BIC rejects real `Position` and neighbouring-Q effects.

A head is the fit module itself (`fit.quality` or `fit.error`): anything with `fit(table, tokens, alphabet)`
and `log_likelihood(components, alphabet, table)`.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from math import log
from typing import Protocol

from sequencing_error_model.observations import CountTable
from sequencing_error_model.spec import Component

CRITERIA = ("aic", "bic", "test-ll")


class Head(Protocol):
    def fit(self, table: CountTable, tokens: Sequence[str], alphabet: Sequence[int]) -> tuple[Component, ...]: ...

    def log_likelihood(self, components: Sequence[Component], alphabet: Sequence[int], table: CountTable) -> float: ...


@dataclass(frozen=True)
class Step:
    step: int
    tokens: tuple[str, ...]
    parameters: int
    test_log_likelihood: float
    criterion: float
    accepted: bool = False


@dataclass(frozen=True)
class Selection:
    components: tuple[Component, ...]
    trace: tuple[Step, ...]


def n_parameters(components: Sequence[Component]) -> int:
    # ponytail: raw array sizes, so softmax gauge and masked cells are overcounted; fine for ranking nested models.
    return sum(v.size for c in components for k, v in c.params.items() if k != "knots")


def select(
    head: Head,
    train: CountTable,
    test: CountTable,
    alphabet: Sequence[int],
    first: Sequence[str],
    candidates: Sequence[Sequence[str]],
    *,
    criterion: str = "test-ll",
) -> Selection:
    """Select head components; `first` lists variants of the required component, each candidate is a group
    of alternative variants (e.g. `("Context(1,1)", "Context(2,2)")`) of which at most one is added."""
    if criterion not in CRITERIA:
        raise ValueError(f"criterion must be one of {CRITERIA}, got {criterion!r}")
    if not first:
        raise ValueError("first must list at least one variant of the required component")
    n = sum(test.counts.values())
    trace: list[Step] = []

    def score(step: int, tokens: tuple[str, ...]) -> tuple[float, tuple[Component, ...]]:
        fitted = head.fit(train, tokens, alphabet)
        ll, k = head.log_likelihood(fitted, alphabet, test), n_parameters(fitted)
        value = {"aic": 2 * k - 2 * ll, "bic": k * log(n) - 2 * ll, "test-ll": -ll}[criterion]
        trace.append(Step(step, tokens, k, ll, value))
        return value, fitted

    def accept(tokens: tuple[str, ...]) -> None:
        i = max(i for i, s in enumerate(trace) if s.tokens == tokens)
        s = trace[i]
        trace[i] = Step(s.step, s.tokens, s.parameters, s.test_log_likelihood, s.criterion, True)

    best_tokens: tuple[str, ...]
    best, best_fit, best_tokens = min(((*score(0, (t,)), (t,)) for t in first), key=lambda x: x[0])
    accept(best_tokens)
    remaining, step = [tuple(c) for c in candidates], 1
    while remaining:
        round_best = None
        for group in remaining:
            for variant in group:
                tokens = (*best_tokens, variant)
                value, fitted = score(step, tokens)
                if round_best is None or value < round_best[0]:
                    round_best = (value, tokens, fitted, group)
        if round_best is None or round_best[0] >= best:
            break
        best, best_tokens, best_fit, chosen = round_best
        remaining.remove(chosen)
        accept(best_tokens)
        step += 1
    return Selection(best_fit, tuple(trace))
