# Analysis plan (pre-registered testing grid)

Registered before the real-data runs in `reports/` (issues #1, #3, #5).

## Testing grid

- **Traits** (2): Alzheimer disease (FinnGen R11 G6_ALZHEIMER; Bellenguez
  2022 GCST90027158 is scripted but could not be staged in the sandbox, see
  docs/DATA_SOURCES.md) and Parkinson disease (FinnGen R11 G6_PARKINSON;
  Nalls 2019 GCST009325 scripted likewise).
- **Hierarchy levels** (3): class (Neuronal/Non-neuronal), supercluster
  (31 WHB superclusters), cluster (461 WHB clusters).
- **Marker rule**: top 100 MapMyCells query markers per taxonomy node
  (protein-coding only; `atlas.build_marker_programs(top_n=100)`).
- **Methods**:
  - MAGMA-style in-repo pipeline (SNP -> gene mean chi-square, 10 kb flank;
    OLS enrichment with log gene SNP-count covariate). **Run in this
    campaign.**
  - Official S-LDSC + MAGMA on the same annotations/programs:
    `scripts/run_ldsc.sh` (1000 Genomes EUR reference required; not
    runnable in the analysis sandbox). **Pending external run.**
- **Multiple testing**: Benjamini-Hochberg FDR at alpha = 0.05 within each
  (trait x level) family.
- **Resolution gain**: number and log10 fold-difference of FDR-significant
  programs at cluster vs supercluster level.
- **Robustness**: sensitivity of supercluster-level z-scores to top_n =
  100 vs 200 markers (Spearman correlation), reported as concordance in lieu
  of cross-method concordance until the S-LDSC run completes.

## Boundary

Statistical association only; no mechanistic claims. Raw GWAS files and raw
atlas expression data are not redistributed; only loaders, checksums, and
small derived tables are committed.
