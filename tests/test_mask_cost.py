from sequencing_error_model import recovery


def test_mask_cost_on_clonal_reads() -> None:
    cost = recovery.mask_cost(recovery.example_spec(), alt_freqs=(0.2,), genome_length=3000, depth=20)
    unmasked, masked = cost["settings"]
    assert unmasked["masked_sites"] == unmasked["removed_rows"] == 0
    assert masked["masked_sites"] > 0 and masked["removed_errors"] > 0
    # Masks select on the outcome: the rows they remove are far more error-rich than the rest.
    assert masked["removed_row_error_rate"] > 2 * masked["row_error_rate"]
