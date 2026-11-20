# Verification log

Independent verification records for the committed Paper 1 results package.
Each entry records only what was actually executed and measured; nothing is
simulated or estimated. Newest entry first.

## 2026-11-20 — independent figure-build and table-integrity verification

**Scope.** Partial re-verification of the results package from committed
inputs only (all of `reports/` as committed at `main`). This covers the
figure leg of #8 and the real-data leg of #7. It is **not** the full
clean-room re-run required by #8: regenerating the enrichment tables from
the frozen input manifest still requires the md5-pinned ABC Atlas download
(~1 GB) and the FinnGen R11 summary statistics, which were not staged in
this environment (same external-data limitation as #10).

**Environment.** Python 3.12.12, matplotlib 3.10.3, pandas 2.3.2, numpy 2.2.5; fresh install of the package with the `dev` extra.

**Checks performed and measured outcomes.**

1. Full test suite: `pytest` — **24 passed** (includes the synthetic-data
   S-LDSC/MAGMA planted-signal recovery tests and the figure determinism
   tests).
2. Real-data figure build:
   `python -m celltype_heritability.figures --reports-dir reports --out-dir figures`
   run twice into independent directories. **11 PNGs written** (volcano +
   top-programs for the 5 committed trait × level cells, plus the
   resolution-gain figure); the PD × cluster cell was skipped with the
   designed warning because its table is not committed (#10). All 11 PNGs
   were **byte-identical across the two runs** (SHA-256 compared).
3. Table integrity vs `reports/RUN_METADATA.json`: every committed
   `magma_enrichment_*` table matches the pinned metadata exactly —
   AD: class 2 tested / 0 FDR, supercluster 30 / 4 (top Oligodendrocyte
   precursor, p = 1.0e-08), cluster 436 / 9 (top `Mgl_9`, p = 2.2e-21);
   PD: class 2 / 0, supercluster 30 / 0 (top Microglia, p = 3.5e-03).
   `resolution_gain.csv` matches the same counts. The PD cluster-level
   table (436 tested, 11 FDR hits, top `CA4_197`) remains uncommitted and
   is tracked in #10.

**Figure digests (this environment).** Run-to-run determinism is the
guaranteed property; absolute bytes can vary with matplotlib/freetype
versions, so these digests identify this environment's build:

| Figure | SHA-256 (first 16 hex) |
|---|---|
| `resolution_gain.png` | `ebe74bd7904e7e27…` |
| `top_programs_AD_FinnGen_R11_class.png` | `38a4dccb75357659…` |
| `top_programs_AD_FinnGen_R11_cluster.png` | `7c5dd71784dc5a0b…` |
| `top_programs_AD_FinnGen_R11_supercluster.png` | `571928b773fa3607…` |
| `top_programs_PD_FinnGen_R11_class.png` | `f0a986b2e612c06d…` |
| `top_programs_PD_FinnGen_R11_supercluster.png` | `9dd9726060723c3f…` |
| `volcano_AD_FinnGen_R11_class.png` | `dfbed85484f46b07…` |
| `volcano_AD_FinnGen_R11_cluster.png` | `14495ac0e964d831…` |
| `volcano_AD_FinnGen_R11_supercluster.png` | `f93778fe1c1fee52…` |
| `volcano_PD_FinnGen_R11_class.png` | `291f950b2b72c407…` |
| `volcano_PD_FinnGen_R11_supercluster.png` | `44043b3873f40fec…` |

**Generated PNGs are not committed**: the automation's GitHub transport
handles text only (binary pushes are corrupted in transit; verified
2026-11-20). Regenerate locally with the command in check 2.
