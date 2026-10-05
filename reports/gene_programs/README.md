# WHB marker gene programs (issue #5)

Deterministic marker programs built from the Allen Brain Cell Atlas whole
human brain (WHB; Siletti et al. 2023) MapMyCells query-marker table by
`celltype_heritability.atlas.build_marker_programs` (top 100 ranked markers
per taxonomy node, protein-coding only, Ensembl -> symbol via the WHB-10Xv3
gene table). Source file versions/md5s are embedded in each JSON under
`source_md5` and documented in `docs/DATA_SOURCES.md`.

## Committed here

* `programs_class_top100.json` - Neuronal / Non-neuronal (2 programs)
* `programs_supercluster_top100.json` - 31 WHB superclusters
* `*_manifest.csv` - one row per program (name, gene count, level)

## Not committed

* `programs_cluster_top100.json` - the 461 cluster-level programs are **not**
  in this directory, but the cluster-level enrichment tables derived from them
  are committed as `reports/magma_enrichment_<trait>_cluster.csv`. Those
  results therefore cannot be regenerated or re-verified from this repository
  alone; the input has to be rebuilt from the Atlas tables below. This is a
  known reproducibility gap, tracked as follow-up work, and no substitute
  file has been written by hand.

Counts are consistent with the committed results: of 31 superclusters, 30 are
testable (the `Bergmann glia` program resolves to 0 mapped protein-coding
markers and is dropped by the `>= 10` gene filter), and 436 of 461 clusters
are testable. See `reports/RUN_METADATA.json`.

Regenerate exactly with:

    python -m celltype_heritability.atlas download --data-dir data/atlas
    python -m celltype_heritability.atlas build --data-dir data/atlas \
        --out-dir reports/gene_programs --level cluster --top-n 100

The `download` step fetches the ABC Atlas tables (~1 GB); the full
MapMyCells precomputed statistics (~8.7 GB) are only needed with `--full`.
