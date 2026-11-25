# Verification log

Independent verification records for the committed Paper 1 results package.
Each entry records only what was actually executed and measured; nothing is
simulated or estimated. Newest entry first.

## 2026-11-25 — publication SVG figure package and committed-table integrity

**Scope.** Automation-reachable follow-up to #7 and the committed-input leg
of #8: extend the figure script from diagnostic panels to the manuscript
figure inventory, commit the text SVG outputs, and move the previously
manual table/metadata checks into the test suite. This is still **not** the
full clean-room re-run required by #8: the external ABC Atlas and FinnGen
inputs remain unstaged in this environment, and the PD cluster-level table
remains the known gap tracked by #10.

**Environment.** Python 3.12.12, matplotlib 3.10.3, pandas 2.3.2,
numpy 2.2.5, scipy 1.16.2; package imported from the local `src/` tree.

**Checks performed and measured outcomes.**

1. Full test suite: `PYTHONPATH=src python -m pytest -q` — **29 passed**
   (24 pre-existing tests plus 5 new figure/integrity checks).
2. Syntax check: `python -m compileall -q src tests scripts` — passed.
3. Real-data figure build:
   `python -m celltype_heritability.figures --reports-dir reports --out-dir <tmp>`
   run twice into independent directories. **28 files written per run**
   (14 unique figures x PNG + SVG): the manuscript Figures 1–3, the
   backwards-compatible resolution-gain panel, and diagnostic volcano /
   top-program panels for the five committed trait x level cells. The PD x
   cluster cell was skipped with the designed warning because its table is
   not committed (#10). All 28 outputs were **byte-identical across the two
   runs** (0 mismatches). Publication SVGs are written by a compact
   fixed-salt-free vector writer; Matplotlib SVG determinism for diagnostic
   panels is pinned by a fixed `svg.hashsalt` and by removing SVG date
   metadata.
4. Committed-table integrity is now CI-covered by
   `tests/test_results_integrity.py`: every committed enrichment table must
   match `RUN_METADATA.json` (row count, FDR count, top program, sorted
   p-values); `resolution_gain.csv` must match the committed tables, with the
   PD cluster cell explicitly sourced from pinned metadata while #10 is open;
   the committed publication SVG paths must exist and be captioned in
   `docs/RESULTS.md`.

**Committed manuscript SVG digests (this environment).** Absolute bytes can
vary with matplotlib/freetype versions, but within this environment the build
is byte-deterministic:

| Figure | SHA-256 (first 16 hex) |
|---|---|
| `figure1_atlas_overview.svg` | `324f81fe397247e0…` |
| `figure2_resolution_gain.svg` | `d8dc8742a3e28a71…` |
| `figure3_AD_FinnGen_R11_cluster_forest.svg` | `0c67cab017137431…` |

No new scientific results were computed; the figures visualize the committed
measured tables only.

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
