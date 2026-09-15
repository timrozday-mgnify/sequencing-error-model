# CLAUDE.md

Guidance for AI agents working in this repo.

- The plan is `docs/IMPLEMENTATION_PLAN.md`. Work in its phase order and update its status markers when a phase lands.
- Tooling is uv + ruff + mypy (strict) + pytest, and CI runs exactly `task lint` and `task test`. Run both before proposing a PR.
- `main` only accepts PRs with green `lint` and `test` checks. Never push to `main` directly.
- Prior art lives in `~/Documents/skiver` (fork of GZHoffie/skiver): `scripts/lib/context_error_models.py`, `error_application.py`, `model_selection.py` and `docs/hmm_error_model.md`. Port ideas and tests from there; don't copy its torch-pickle artifact format.
- **Default mode must work with an unmodified, pinned skiver release.** Never make a default-mode code path depend on fork-only outputs (`skiver dump`, `windows.bin`, `fragment_id`, `read_id`).
- Every exporter needs a round-trip test: export, simulate, re-estimate, compare rates within tolerance.
- Don't commit large skiver outputs or trained models. Test fixtures must be small (under 1 MB; the pre-commit hook enforces this).
