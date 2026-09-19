# GWAS × Brain Cell-Type Heritability (Paper 1)

This independent research repository plans and tracks a cell-type-resolved heritability analysis of Alzheimer's and Parkinson's disease GWAS across the modern whole-brain cell-type hierarchies. It provides data-free analysis utilities for transparent review and extension.

## Series position

This is **Paper 1** of the GWAS × brain cell-type series (4 papers). It is the foundation of the series; Papers 2–4 build on its outputs.

## Research plan

| Planned work | Expected outcome |
|---|---|
| Assemble open GWAS summary statistics (AD: Bellenguez 2022; PD: Nalls 2019/Kim 2024) | Versioned, harmonized summary-statistics inputs. |
| Build cell-type gene programs from Allen ABC Atlas human whole brain (~3M nuclei, >3,000 types) | Cell-type/supertype marker sets at multiple hierarchy levels. |
| Run S-LDSC / MAGMA / sc-linker across hierarchy levels | Heritability enrichment atlas per trait × cell-type level. |
| Compare fine vs. coarse resolutions | Quantified gain of the new fine-grained taxonomies over legacy ~8-type labels. |

**Current status:** statistics plumbing implemented and validated on synthetic data; all inputs are open data; no real-data analysis has been run.

## What is included

| Path | Contents |
|---|---|
| `src/celltype_heritability/` | Data-free enrichment utilities: `annotations.py` (binary cell-type SNP annotations), `ldscore.py` (toy partitioned LD scores), `sldsc.py` (S-LDSC-style WLS with jackknife SEs), `magma.py` (gene-level association + program enrichment), `sclinker.py` (disease-program scoring skeleton), `simulate.py` (synthetic GWAS with planted enrichment). |
| `tests/` | Synthetic tests: planted cell-type enrichment is recovered by both the S-LDSC and MAGMA paths (`python -m pytest -q`). |
| `docs/` | Research status, methods scope, contribution guidance, and [data-access prerequisites](docs/DATA_ACCESS.md). |
| `.github/` | CI workflow (pytest, Python 3.10–3.12), task-card issue template, PR template. |

## Use and validation

```bash
pip install -e ".[dev]"
python -m pytest -q
```

## Keywords

GWAS, Alzheimer's disease, Parkinson's disease, S-LDSC, MAGMA, sc-linker, single-cell, Allen Brain Atlas, neurogenetics, reproducible research.

## Documentation

- [Introduction for new readers](docs/INTRODUCTION.md)
- [Data access and prerequisites](docs/DATA_ACCESS.md)
- [Contributing](CONTRIBUTING.md)
