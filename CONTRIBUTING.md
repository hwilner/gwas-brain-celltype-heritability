# Contributing

Thanks for your interest in the GWAS × brain cell-type heritability project
(Paper 1 of 4). Contributions of all sizes are welcome.

## Workflow

1. Pick a task card from the issue tracker. Cards carry `task-key` metadata,
   size labels (`S`/`M`/`L`), and acceptance criteria — read them fully
   before starting.
2. Open a pull request against `main` using the PR template. Reference the
   issue number and restate the acceptance criteria you are meeting.
3. CI (pytest on Python 3.10–3.12) must be green before review.

## Development setup

```bash
pip install -e ".[dev]"
python -m pytest -q
```

All statistics plumbing must be validated on **synthetic data** (see
`src/celltype_heritability/simulate.py`); tests must not require restricted
or large external datasets.

## Scope rules

- Code and synthetic-data tests can be merged once CI is green.
- Tasks gated by `gate: scientific review` (e.g. real-data runs, manuscript)
  additionally require scientific-owner sign-off.
- Do not commit restricted data or redistribute controlled-access files;
  record sources and checksums in `docs/` instead.
