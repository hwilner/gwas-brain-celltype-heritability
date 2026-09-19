# Introduction — Paper 1: Cell-Type-Resolved Heritability Atlas of AD and PD

**Series note:** This is **Paper 1 of 4** in the GWAS × brain cell-type series. It is the first paper and does not build on any former paper. Papers 2 (cross-species mapping), 3 (coloc/TWAS prioritization), and 4 (psychiatric extension) all reuse the harmonized GWAS inputs and cell-type gene programs built here.

## Background

Cell-type enrichment analyses of brain GWAS have so far used coarse (~8) cortical cell-type labels. The 2023–24 whole-brain atlases (Allen ABC Atlas human: ~3M nuclei, >3,000 types; mouse: 1,201 supertypes) allow heritability to be localized at far finer resolution — including glial *states*, not just types. All required inputs are fully open: GWAS summary statistics (EBI GWAS Catalog, PGC) and processed atlas data (ABC Atlas, CELLxGENE).

## Research questions

1. Which human brain cell types/supertypes mediate AD and PD GWAS heritability?
2. Do AD loci converge on microglial/astrocyte states rather than types?
3. How much resolution is gained over legacy coarse labels?

## Data

| Resource | Scale | Access |
|---|---|---|
| AD GWAS (Bellenguez 2022, GCST90027158) | ~789k samples | Open (GWAS Catalog) |
| PD GWAS (Nalls 2019 / Kim 2024) | large meta-analyses | Open |
| Allen ABC Atlas human whole brain | ~3M nuclei | Open (AWS) |
| SEA-AD disease-state annotations | 84 donors | Open |

## Methods

S-LDSC stratified heritability, MAGMA gene-set analysis, sc-linker enhancer-gene linking; multiple hierarchy levels; FDR control.

## Expected contributions

- A fine-grained heritability atlas of AD/PD across the brain cell-type hierarchy.
- Versioned, reusable gene-program resources for Papers 2–4.

## Scope and boundary

Planning, software, and synthetic tests live here. Empirical outputs are produced under data-use compliance and released upon owner decision. Enrichment is statistical association, not mechanism.
