# WHB marker gene programs (issue #5)

Deterministic marker programs built from the Allen Brain Cell Atlas whole
human brain (WHB; Siletti et al. 2023) MapMyCells query-marker table by
`celltype_heritability.atlas.build_marker_programs` (top 100 ranked markers
per taxonomy node, protein-coding only, Ensembl -> symbol via the WHB-10Xv3
gene table). Source file versions/md5s are embedded in each JSON under
`source_md5` and documented in `docs/DATA_SOURCES.md`.

* `programs_class_top100.json` - Neuronal / Non-neuronal (2 programs)
* `programs_supercluster_top100.json` - 31 WHB superclusters
* `programs_cluster_top100.json` - 461 WHB clusters (finest level)
* `*_manifest.csv` - one row per program (name, gene count, level)

Regenerate exactly with:

    python -m celltype_heritability.atlas download --data-dir data/atlas
    python -m celltype_heritability.atlas build --data-dir data/atlas \
        --out-dir reports/gene_programs --level cluster --top-n 100
