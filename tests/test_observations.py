from collections import Counter

import pytest

from sequencing_error_model.observations import CountTable, subset


def test_marginal() -> None:
    t = CountTable("x", ("q", "op"), "base", True, Counter({(2, "="): 5, (2, "!"): 1, (37, "="): 9}))
    assert t.marginal("op").counts == {("=",): 14, ("!",): 1}
    assert t.marginal().counts == {(): 15}


def test_expected_counts_are_fractional() -> None:
    t = CountTable("x", ("q", "op"), "base", True, Counter({(2, "="): 0.25, (2, "A>C"): 0.75}))
    assert t.marginal("q").counts == {(2,): 1.0}


def test_op_needs_truth() -> None:
    with pytest.raises(ValueError, match="quality is never error truth"):
        CountTable("x", ("q", "op"), "base", False, Counter({(2, "!"): 1}))


@pytest.mark.parametrize(
    ("fields", "counts", "match"),
    [
        (("q", "qual"), {(2, 3): 1}, "unknown or repeated"),
        (("q", "q"), {(2, 3): 1}, "unknown or repeated"),
        (("q", "op"), {(2,): 1}, "bad entry"),
        (("q", "op"), {(2, "="): -1}, "bad entry"),
        (("q", "op"), {(2, "A>A"): 1}, "unknown op"),
    ],
)
def test_invalid(fields: tuple[str, ...], counts: dict[tuple[object, ...], int], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        CountTable("x", fields, "base", True, Counter(counts))


def test_subset_strides() -> None:
    assert list(subset(range(10))) == list(range(10))
    assert list(subset(range(10), every_n=3)) == [0, 3, 6, 9]
    assert list(subset(range(10), 2, 3)) == [0, 3]
    assert list(subset(range(10), every_n=3, offset=1)) == [1, 4, 7]
    # Disjoint samples, so a second one can check the first.
    assert not set(subset(range(10), every_n=3)) & set(subset(range(10), every_n=3, offset=1))


@pytest.mark.parametrize(("every_n", "offset"), [(0, 0), (-1, 0), (1, -1)])
def test_subset_invalid(every_n: int, offset: int) -> None:
    with pytest.raises(ValueError, match="every_n must be"):
        subset(range(3), every_n=every_n, offset=offset)
