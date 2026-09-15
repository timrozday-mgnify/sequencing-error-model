# Contributing

## Setup

```bash
uv sync
uv run pre-commit install --hook-type pre-commit --hook-type pre-push
```

## PR policy

`main` is protected by a GitHub ruleset:

- No direct pushes: every change lands through a pull request.
- Required status checks, which must pass on the latest commit: `lint` (Linting workflow) and `test` (Testing workflow).
- No force-pushes to or deletion of `main`.

Branch names are short and descriptive (`skiver-csv-parser`, `export-badread`).
Keep PRs to one plan phase or one feature. Fill in the PR template. It asks for
an explicit note on any change to an on-disk format (model artifact or exporter
output), because downstream simulators consume those files.

## Local checks (same as CI)

```bash
task lint   # ruff format --check, ruff check, mypy
task test   # pytest
```

## Releases

Bump `version` in `pyproject.toml`, merge, then tag `vX.Y.Z` on `main` and push the tag.
The Release workflow checks that the tag matches the package version, runs the tests,
builds, and publishes a GitHub release.
