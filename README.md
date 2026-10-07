# GWAS × Brain Cell-Type Heritability 

This independent research repository plans and tracks a cell-type-resolved heritability analysis of Alzheimer's and Parkinson's disease GWAS across the modern whole-brain cell-type hierarchies. It provides data-free analysis utilities for transparent review and extension.

## Introduction

Common-variant heritability of Alzheimer (AD) and Parkinson (PD) disease is concentrated in regulatory DNA, and single-cell genomics now offers whole-brain taxonomies at unprecedented resolution (Allen Brain Cell Atlas whole human brain: ~3M nuclei, >3,000 types; Siletti et al. 2023). This project tests which cell types carry GWAS signal for each trait, and whether the new fine-grained taxonomies recover signal that legacy coarse labels miss. See [docs/INTRODUCTION.md](docs/INTRODUCTION.md) for the full introduction for new readers.

## Extended introduction

Background, prior art (S-LDSC, MAGMA, sc-linker), the ABC Atlas taxonomy, and the rationale for the pre-registered design are developed in [docs/EXTENDED_INTRODUCTION.md](docs/EXTENDED_INTRODUCTION.md). Data sources and checksums are pinned in [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md); access prerequisites in [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md).

## Materials and methods

Full methods, including the explicit done-vs-intended separation, are in [docs/METHODS.md](docs/METHODS.md); the pre-registered testing grid (2 traits × 3 hierarchy levels, top-100 markers, BH-FDR α = 0.05 per family) is in [docs/ANALYSIS_PLAN.md](docs/ANALYSIS_PLAN.md). In brief:

| Component | Implementation |
|---|---|
| Cell-type gene programs | ABC Atlas WHB MapMyCells top-100 ranked markers per node, protein-coding only, at class / supercluster / cluster levels, md5-pinned (`reports/gene_programs/`, `atlas.py`) |
| GWAS inputs | FinnGen R11 `G6_ALZHEIMER` and `G6_PARKINSON`, deterministic QC (21.3M → 18.57M variants), md5-pinned (`sumstats.py`) |
| Enrichment tests | In-repo MAGMA-style pipeline (SNP→gene ±10 kb, mean chi-square, OLS program enrichment with log-SNP-count covariate; `magma.py`) — run; official S-LDSC scripted (`scripts/run_ldsc.sh`) — pending external LD reference |
| Validation | Synthetic GWAS with planted enrichment, recovered by both S-LDSC and MAGMA paths (`tests/`) |

## Results

Captioned paper-style tables for all committed runs are in [docs/RESULTS.md](docs/RESULTS.md); machine-readable tables and pinned run metadata are in `reports/`. Headline findings:

**Table 1. Enrichment atlas headline results (in-repo MAGMA-style pipeline, FDR < 0.05 per trait × level family).** Fine cluster resolution recovers signal invisible at coarse levels — most sharply for PD, which has zero supercluster-level but eleven cluster-level hits.

| Trait | class (2 tested) | supercluster (30 tested) | cluster (436 tested) | Top hit |
|---|---|---|---|---|
| AD | 0 | 4 (OPC q = 3.1e-07; Vascular; Microglia; Astrocyte) | 9 (glia-dominated) | `Mgl_9`, p = 2.2e-21 |
| PD | 0 | 0 (top: Microglia, q = 0.106) | 11 | `CA4_197`, p = 5.3e-13 |

**Current status:** in-repo MAGMA-style runs on FinnGen R11 AD/PD are complete and committed (Tables above; robustness: top-100 vs top-200 marker Spearman ρ = 0.81/0.86). Official S-LDSC and EBI-staged GWAS replications are scripted but not yet run — see open issues.

## Research plan

| Planned work | Expected outcome |
|---|---|
| Assemble open GWAS summary statistics (AD: Bellenguez 2022; PD: Nalls 2019/Kim 2024) | Versioned, harmonized summary-statistics inputs. |
| Build cell-type gene programs from Allen ABC Atlas human whole brain (~3M nuclei, >3,000 types) | Cell-type/supertype marker sets at multiple hierarchy levels. |
| Run S-LDSC / MAGMA / sc-linker across hierarchy levels | Heritability enrichment atlas per trait × cell-type level. |
| Compare fine vs. coarse resolutions | Quantified gain of the new fine-grained taxonomies over legacy ~8-type labels. |

## What is included

| Path | Contents |
|---|---|
| `src/celltype_heritability/` | Data-free enrichment utilities: `annotations.py` (binary cell-type SNP annotations), `ldscore.py` (toy partitioned LD scores), `sldsc.py` (S-LDSC-style WLS with jackknife SEs), `magma.py` (gene-level association + program enrichment), `sclinker.py` (disease-program scoring skeleton), `simulate.py` (synthetic GWAS with planted enrichment). |
| `tests/` | Synthetic tests: planted cell-type enrichment is recovered by both the S-LDSC and MAGMA paths (`python -m pytest -q`). |
| `docs/` | Research status, methods scope, results, contribution guidance, and [data-access prerequisites](docs/DATA_ACCESS.md). |
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
- [Extended introduction](docs/EXTENDED_INTRODUCTION.md)
- [Methods](docs/METHODS.md) · [Analysis plan](docs/ANALYSIS_PLAN.md) · [Results](docs/RESULTS.md)
- [Data access and prerequisites](docs/DATA_ACCESS.md)
- [Contributing](CONTRIBUTING.md)
