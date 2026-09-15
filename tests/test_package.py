import sequencing_error_model


def test_version_matches_distribution() -> None:
    assert sequencing_error_model.__version__ == "0.0.1"
