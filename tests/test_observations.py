from collections import Counter

import pytest

from sequencing_error_model.observations import CountTable


def test_marginal() -> None:
    t = CountTable("x", ("q", "op"), "base", True, Counter({(2, "="): 5, (2, "!"): 1, (37, "="): 9}))
    assert t.marginal("op").counts == {("=",): 14, ("!",): 1}
    assert t.marginal().counts == {(): 15}


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
