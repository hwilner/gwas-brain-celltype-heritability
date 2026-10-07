# Results

*Paper-style report of the committed real-data runs. All numbers are taken
verbatim from `reports/` (run metadata: `reports/RUN_METADATA.json`);
the pre-registered testing grid is in [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md)
and the methods in [METHODS.md](METHODS.md). Statistical association only;
no mechanistic claims.*

## Run overview

**Table 1. Analysis runs committed in this repository.** In-repo MAGMA-style
pipeline (`scripts/run_magma_enrichment.py`): SNPs mapped to genes with a
±10 kb flank, gene statistic = mean SNP chi-square, per-program OLS
enrichment with a log gene SNP-count covariate; Benjamini–Hochberg FDR at
α = 0.05 within each (trait × level) family. Cell-type programs are the
top-100 ranked MapMyCells markers per Allen WHB taxonomy node
(protein-coding only). The official S-LDSC run on the same programs is
scripted (`scripts/run_ldsc.sh`) but pending an external LD reference, so
no S-LDSC results are reported here.

| Trait | GWAS source | Variants in → out (QC) | Genes tested | Levels tested (programs) | FDR < 0.05 hits |
|---|---|---|---|---|---|
| Alzheimer disease | FinnGen R11 `G6_ALZHEIMER` (11,755 cases / 441,978 controls) | 21,306,794 → 18,565,791 | 18,506 | class (2), supercluster (30), cluster (436) | 0 / 4 / 9 |
| Parkinson disease | FinnGen R11 `G6_PARKINSON` (5,150 cases / 448,583 controls) | 21,306,794 → 18,565,809 | 18,506 | class (2), supercluster (30), cluster (436) | 0 / 0 / 11 |

Input files are md5-pinned in `reports/RUN_METADATA.json` (GWAS and ABC
Atlas tables); of 31 superclusters, 30 are testable (the Bergmann-glia
program resolves to 0 mapped protein-coding markers), and 436 of 461
clusters pass the ≥10-gene filter.

## Alzheimer disease

**Table 2. AD heritability enrichment, supercluster level — all
FDR-significant programs (4 of 30 tested).** Full table:
`reports/magma_enrichment_AD_FinnGen_R11_class.csv` (class level, 0 of 2
significant) and `reports/magma_enrichment_AD_FinnGen_R11_supercluster.csv`.
β is the program-membership coefficient from the enrichment regression;
q is the Benjamini–Hochberg adjusted p-value within the AD × supercluster
family.

| Program | β | SE | z | p | q (FDR) |
|---|---|---|---|---|---|
| Oligodendrocyte precursor | 1.275 | 0.223 | 5.73 | 1.0e-08 | 3.1e-07 |
| Vascular | 1.198 | 0.218 | 5.51 | 3.7e-08 | 5.5e-07 |
| Microglia | 0.838 | 0.218 | 3.84 | 1.2e-04 | 1.2e-03 |
| Astrocyte | 0.771 | 0.219 | 3.52 | 4.4e-04 | 3.3e-03 |

All remaining 26 superclusters are null (minimum q among them = 0.88). The
coarse class level (Neuronal vs Non-neuronal) shows no signal (both
q = 0.91), i.e. the enrichment is invisible at legacy ~2-class resolution.

**Table 3. AD heritability enrichment, cluster level — all FDR-significant
programs (9 of 436 tested).** Full table:
`reports/magma_enrichment_AD_FinnGen_R11_cluster.csv`. Glial programs
dominate: four microglial (`Mgl_*`), one oligodendrocyte-precursor
(`OPC_32`), one pericyte (`Per_21`), and two astrocyte (`Astro_*`) clusters;
`Midi_435` is a midbrain-derived inhibitory neuronal cluster.

| Program | β | SE | z | p | q (FDR) |
|---|---|---|---|---|---|
| `Mgl_9` | 5.276 | 0.556 | 9.50 | 2.2e-21 | 9.5e-19 |
| `Mgl_12` | 2.536 | 0.291 | 8.72 | 2.8e-18 | 6.2e-16 |
| `OPC_32` | 2.948 | 0.415 | 7.10 | 1.3e-12 | 1.8e-10 |
| `Per_21` | 1.506 | 0.223 | 6.74 | 1.5e-11 | 1.4e-09 |
| `Midi_435` | 1.695 | 0.252 | 6.73 | 1.7e-11 | 1.4e-09 |
| `Mgl_6` | 1.295 | 0.226 | 5.74 | 9.6e-09 | 7.0e-07 |
| `Mgl_7` | 1.882 | 0.349 | 5.39 | 7.2e-08 | 4.5e-06 |
| `Astro_53` | 1.287 | 0.289 | 4.45 | 8.7e-06 | 4.7e-04 |
| `Astro_60` | 0.764 | 0.219 | 3.49 | 4.8e-04 | 2.3e-02 |

## Parkinson disease

**Table 4. PD heritability enrichment summary.** No program passes FDR at
class (0 of 2; both q = 0.54) or supercluster (0 of 30) level; the top
supercluster hit, Microglia, is nominally significant only. At cluster
level 11 of 436 programs pass FDR < 0.05 (top: `CA4_197`, minimum
p = 5.3e-13). Source: `reports/magma_enrichment_PD_FinnGen_R11_class.csv`,
`reports/magma_enrichment_PD_FinnGen_R11_supercluster.csv`,
`reports/RUN_METADATA.json`, `reports/resolution_gain.csv`.

| Level | Programs tested | FDR < 0.05 | Top program | Top-program statistic |
|---|---|---|---|---|
| class | 2 | 0 | Neuronal | β = −0.032, p = 0.50 |
| supercluster | 30 | 0 | Microglia | β = 0.132, z = 2.92, p = 3.5e-03, q = 0.106 |
| cluster | 436 | 11 | `CA4_197` | p = 5.3e-13 (minimum p) |

> **Known gap:** the PD cluster-level full table
> (`magma_enrichment_PD_FinnGen_R11_cluster.csv`) is not committed; the
> counts above come from the pinned run metadata. Regenerating it requires
> rebuilding the cluster-level programs from the md5-pinned Atlas tables
> (see `reports/gene_programs/README.md`). Tracked as issue #10.

## Resolution gain

**Table 5. Gain of fine (cluster) over coarse (supercluster) taxonomy
resolution.** Counts of FDR-significant programs and minimum p-value per
level. Source: `reports/resolution_gain.csv`.

| Trait | Sig. supercluster | Sig. cluster | Min p (supercluster) | Min p (cluster) |
|---|---|---|---|---|
| AD | 4 | 9 | 1.0e-08 | 2.2e-21 |
| PD | 0 | 11 | 3.5e-03 | 5.3e-13 |

For PD the entire signal is a resolution effect: zero supercluster-level
hits versus eleven cluster-level hits.

## Robustness

**Table 6. Sensitivity of program z-scores to the marker-count rule.**
Spearman correlation between supercluster-level z-scores computed with
top-100 vs top-200 markers per node. Source:
`reports/concordance_topn_sensitivity.csv`. This stands in for cross-method
concordance until the official S-LDSC run completes (issue #8).

| Trait | Spearman ρ (top-100 vs top-200 z) |
|---|---|
| AD | 0.81 |
| PD | 0.86 |

## Figures

No results figures are committed yet; figure scripts are tracked as issue
#7. Concept schematics (`01-concept-schematic.png`,
`02-data-collection.png`) live in `docs/figures/`.

## Verification record

The synthetic-data test suite (planted enrichment recovered by both the
S-LDSC and MAGMA paths) passes on all branches — see the audit table in the
programme state file and `tests/`. Unit tests verify code behaviour only;
they do not verify the scientific claims above, which rest on the pinned
inputs and committed tables in `reports/`.
