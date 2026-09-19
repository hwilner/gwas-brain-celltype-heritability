# Data access

> **Status: placeholder.** This page lists the external inputs required
> before the real-data analyses (issues #1, #3, #5) can run. Nothing here is
> redistributed; all inputs are downloaded by the user from public sources.

## Allen Brain Cell (ABC) Atlas

- Resource: Allen Brain Cell Atlas, human whole-brain snRNA-seq/snATAC-seq
  (~3M nuclei, >3,000 cell types).
- Download: https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas
  (AWS open-data bucket; see the portal for `abc_atlas_access` instructions).
- Needed for: cell-type marker genes and cell-type-specific peaks used to
  build gene programs and SNP annotations at each hierarchy level
  (class → subclass → supertype/cluster).

## LD reference panel (for S-LDSC)

- 1000 Genomes Phase 3 European reference panel and pre-computed partitioned
  LD scores / regression weights, from the official `ldsc` resources:
  https://github.com/bulik/ldsc (see the wiki; `1000G_EUR_Phase3` files).
- Production S-LDSC runs use the official `ldsc` software; the in-repo
  `ldscore.py` is a toy implementation for synthetic tests only.

## GWAS summary statistics (public)

- Alzheimer's disease: Bellenguez et al. 2022 (EBI GWAS Catalog
  GCST90027158).
- Parkinson's disease: Nalls et al. 2019 (excluding 23andMe) and/or
  Kim et al. 2024 (GP2 / public release).
- Harmonization targets a common build (GRCh38) and format; sources,
  versions, and checksums will be recorded in `docs/DATA_SOURCES.md`
  (issue #3). Restricted files (e.g. 23andMe) are never redistributed.

## Once staged

Place downloaded files under `data/` (git-ignored) and the analysis code in
`src/celltype_heritability/` is ready to run against them; see the
`simulate.py` docstrings for the expected schemas.
