import re
from dataclasses import replace

import numpy as np

from sequencing_error_model import compare, recovery
from sequencing_error_model import generate as gen
from sequencing_error_model.fit import error
from sequencing_error_model.sources import pe_overlap

PCR = 0.01  # library substitutions per fragment base, shared by both mates


def test_pe_overlap_vs_reference_estimates_library_errors() -> None:
    example, rng = recovery.example_spec(), np.random.default_rng(7)
    head0 = example.error_head[0]
    bias = head0.params["bias"] - np.r_[0, np.full(9, 1.5)]  # ~3.5% errors, as in the pe-overlap recovery test
    truth = replace(
        example, error_head=(replace(head0, params={**head0.params, "bias": bias}), *example.error_head[1:])
    )
    genome, length = "".join(rng.choice(list("ACGT"), size=20_000)), 40

    def ends(frag: str) -> list[str]:
        return [
            (f + a + "N" * length)[:length] for f, a in ((frag, gen.ADAPTERS[0]), (gen._revcomp(frag), gen.ADAPTERS[1]))
        ]

    true_templates, read_templates = [], []
    for name, t1, t2 in gen.fragments([("g", genome)], gen.insert_sizes(45, 8), 3000, length, rng):
        match = re.fullmatch(r"g:(\d+)-(\d+)#\d+", name)
        assert match
        frag = genome[int(match[1]) - 1 : int(match[2])]
        assert ends(frag) == [t1, t2]
        hit = rng.random(len(frag)) < PCR
        library = "".join(
            str(rng.choice([x for x in "ACGT" if x != c])) if h else c for c, h in zip(frag, hit, strict=True)
        )
        true_templates += [t1, t2]
        read_templates += ends(library)
    mates = [1, 2] * (len(true_templates) // 2)
    reads = gen.generate(truth, read_templates, mates, rng)

    # Both mates read the library's substitutions, so they cancel in the overlap and show against the genome.
    pairs = (
        (a.sequence, [ord(c) - 33 for c in a.quality], b.sequence, [ord(c) - 33 for c in b.quality])
        for a, b in zip(reads[::2], reads[1::2], strict=True)
    )
    ev = pe_overlap.collect(pairs, 1, (2, 2))
    tokens = [c.token for c in truth.error_head]
    overlap_head = pe_overlap.fit(ev, tokens, ["QualityMarkov(1)"], iterations=5).error_head
    overlap = pe_overlap.table(ev, overlap_head)
    reference = gen.observations(zip(true_templates, reads, mates, strict=True), (2, 2), 1, "reference")

    found = compare.evidence(overlap, reference, by=("q", "mate"))
    assert abs(found["excess"] - PCR) < 0.2 * PCR, found

    # Model level: the reference-fitted head E carries the library errors on top of the overlap's.
    spec_a = replace(truth, error_head=overlap_head)
    spec_b = replace(truth, error_head=error.fit(reference, tokens, truth.quality_alphabet))
    fitted = compare.models(spec_a, spec_b, reference)
    assert abs(fitted["excess"] - PCR) < 0.3 * PCR, fitted
    assert fitted["log_likelihood_per_row"][1] > fitted["log_likelihood_per_row"][0]
