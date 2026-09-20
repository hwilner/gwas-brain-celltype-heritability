# Data sources

All inputs are public; no restricted data (e.g. 23andMe) is used or
redistributed. Raw files are downloaded to `data/` (git-ignored) by the
loaders in `src/celltype_heritability/`; only small derived tables are
committed under `reports/`.

## Cell-type atlas (issue #5)

Allen Brain Cell Atlas, human whole-brain (WHB) taxonomy of Siletti et al.
2023 (Science, PMID 37676928), 3.3M nuclei, 461 clusters, 31 superclusters.
Marker genes are the published MapMyCells query-marker table
(`cell_type_mapper` v1.3.0, computed 2024-02-21) served from the ABC Atlas
open S3 bucket (release 20241130):

| file | url (prefix https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/) | md5 |
|---|---|---|
| query_markers.n10.20240221800.json | mapmycells/WHB-10Xv3/20240831/query_markers.n10.20240221800.json | 286e1fcc3491b11a13ce7cb35a618733 |
| cluster.csv | metadata/WHB-taxonomy/20240330/cluster.csv | c5648d1c646f64c590ff9c9dac53e17b |
| cluster_annotation_term.csv | metadata/WHB-taxonomy/20240330/cluster_annotation_term.csv | 2533ac1498b35a1ce386a3d97137e9c3 |
| cluster_annotation_term_set.csv | metadata/WHB-taxonomy/20240330/cluster_annotation_term_set.csv | 89d131f4afabf925019e4d405c2a150d |
| cluster_to_cluster_annotation_membership.csv | metadata/WHB-taxonomy/20240330/cluster_to_cluster_annotation_membership.csv | ae0850e319e8a3162d72bdcc971106e1 |
| gene.csv (WHB-10Xv3) | metadata/WHB-10Xv3/20241115/gene.csv | d3400fda8166398e49e994cf36a39028 |

Optional full reference (`--full` flag of
`python -m celltype_heritability.atlas download`): MapMyCells
`precomputed_stats.siletti.training.h5`, ~8.7 GB (streamed in chunks;
not needed for the committed analyses).

Loader: `celltype_heritability.atlas.download_atlas` /
`build_marker_programs`. Programs committed in `reports/gene_programs/`
(top 100 ranked markers per node, protein-coding only; class / supercluster /
cluster levels, with manifests).

## Gene annotation

Ensembl GRCh38 gene table via BioMart (`useast.ensembl.org/biomart`), gene
stable ID / symbol / chromosome / span / biotype; protein-coding autosomal
genes, longest span per symbol -> 18,551 genes
(`sumstats.load_gene_annotation`).

## GWAS summary statistics (issue #3)

### Staged and analyzed (real)

| trait | source | variants in file | after QC | md5 |
|---|---|---|---|---|
| AD | FinnGen R11 `finngen_R11_G6_ALZHEIMER.gz` (GRCh38; 11,755 cases / 441,978 controls) | 21,306,794 | 18,565,791 | 30a0a6701918fedfe84665cf33e46eb2 |
| PD | FinnGen R11 `finngen_R11_G6_PARKINSON.gz` (GRCh38; 5,150 cases / 448,583 controls) | 21,306,794 | 18,565,809 | da8da233c84966b069b8a9f2fd154363 |

URLs: `https://storage.googleapis.com/finngen-public-data-r11/summary_stats/`
(public FinnGen release 11 bucket). QC (deterministic,
`sumstats.stream_harmonize`): autosomes 1-22, 0 < p <= 1, non-missing
beta/se with se > 0, ACGT SNP alleles, deduplicated by (chrom, pos).

### Documented and scripted, not staged in this sandbox

The issue-#3 named datasets require downloading 0.5-0.8 GB from the EBI GWAS
Catalog FTP, which was throttled to ~40 KB/s in this environment (multi-hour
transfers). They are registered in `sumstats.GWAS_SOURCES` with working
resumable downloaders; official checksums from the EBI `md5sum.txt` files:

| trait | accession | file | size | official md5 |
|---|---|---|---|---|
| AD | GCST90027158 (Bellenguez 2022) | GCST90027158_buildGRCh38.tsv.gz | 755,201,909 B | 9d23b9ba23532da38ab83fb061bab18f |
| PD | GCST009325 (Nalls 2019 excl. 23andMe) | harmonised/GCST009325.h.tsv.gz | 490,619,051 B | e7f2654a34931bb9996fb0927cece282 |

Fetch with: `python -m celltype_heritability.sumstats download
AD_Bellenguez2022 --dest data/gwas/` (resumable). Re-running
`scripts/run_magma_enrichment.py` with these files reproduces the atlas
against the larger GWAS.

## LD reference (for official S-LDSC)

1000 Genomes Phase 3 EUR (`1000G_Phase3_plinkfiles`,
`1000G_Phase3_baselineLD_v2.2_ldscores`, `weights_hm3_no_hla`) from
https://data.broadinstitute.org/alkesgroup/LDSCORE/. Required by
`scripts/run_ldsc.sh` (official `ldsc`/`magma` production run; not executed
in this sandbox - see `docs/ANALYSIS_PLAN.md` for what is real vs pending).
