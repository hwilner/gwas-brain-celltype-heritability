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

**Current status:** planning stage; all inputs are open data; no analysis has been run.

## What is included

| Path | Contents |
|---|---|
| `src/` | In-memory enrichment-analysis utilities (added as issues are completed). |
| `tests/` | Synthetic tests for enrichment statistics. |
| `docs/` | Research status, methods scope, and contribution guidance. |

## Use and validation

```bash
python -m pytest -q
```

## Keywords

GWAS, Alzheimer's disease, Parkinson's disease, S-LDSC, MAGMA, sc-linker, single-cell, Allen Brain Atlas, neurogenetics, reproducible research.

## Documentation

- [Introduction for new readers](docs/INTRODUCTION.md)
