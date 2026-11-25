# Paper 1 — Manuscript Outline

*Working title:* **Fine-grained brain cell-type taxonomies recover GWAS signal invisible to coarse labels: an Alzheimer's and Parkinson's enrichment atlas**

This outline is the reviewer-facing skeleton of Paper 1. It states the paper's argument, its position in the series, and the exact table/figure inventory the finished manuscript must contain. Section drafts draw only on committed, versioned results (`reports/`, `docs/RESULTS.md`); no claim appears here that the repository cannot already back.

---

## Series position

Paper 1 of 4 in the GWAS × brain cell-type series. This paper is the **foundation**: it builds the harmonized inputs, the cell-type gene programs, and the enrichment testing grid that every later paper reuses.

### Resources this paper produces for Papers 2–4

| Resource | Where it lives | Reused by |
|---|---|---|
| Harmonized, md5-pinned FinnGen R11 GWAS inputs with deterministic QC (21.3M → 18.57M variants) | `reports/RUN_METADATA.json`, `sumstats.py` | Papers 2–4 |
| ABC Atlas WHB cell-type gene programs (top-100 MapMyCells markers, protein-coding, class/supercluster/cluster levels) | `reports/gene_programs/`, `atlas.py` | Papers 2–4 |
| In-repo MAGMA-style enrichment pipeline + scripted S-LDSC harness | `magma.py`, `scripts/run_magma_enrichment.py`, `scripts/run_ldsc.sh` | Papers 2–4 |
| Pre-registered testing grid (2 traits × 3 levels, BH-FDR α = 0.05 per family) | `docs/ANALYSIS_PLAN.md` | Papers 2–4 |
| Prioritized cell types and trait associations (this atlas) | `docs/RESULTS.md` | Papers 2–3 (downstream prioritization and cross-species concordance) |

Papers 2–4 depend on this paper's frozen inputs and program tables; a change to the gene-program build after publication is a versioning event, not an edit.

---

## The argument in three sentences

1. Common-variant signal for AD and PD is concentrated in regulatory DNA, but standard cell-type enrichment analyses label cells with legacy coarse taxonomies (~2–8 types).
2. Re-running the analysis against the Allen whole-human-brain hierarchy (3,000+ types) shows the coarse labels *hide* the signal: AD shows nothing at class level but four supercluster and nine cluster hits; **PD shows nothing at class *or* supercluster level, and eleven cluster-level hits** — its entire detectable signal is a resolution effect.
3. Fine-grained taxonomies are therefore not a refinement but a requirement for this class of analysis, and the atlas this paper releases is the substrate for locus-to-gene and cross-species follow-up (Papers 2–4).

---

## Section structure

### 1. Introduction
- GWAS heritability concentrates in regulatory DNA; the cell-type-resolution question.
- Prior art and its ceiling: S-LDSC, MAGMA, sc-linker run against legacy ~8-type labels.
- The ABC Atlas whole-human-brain taxonomy (Siletti et al. 2023; ~3M nuclei, >3,000 types) as the new substrate.
- This paper's pre-registered grid and its discipline: statistical association only, no mechanistic claims; FDR controlled per trait × level family.
- Series framing: Paper 1 of 4; explicitly name the resources Papers 2–4 reuse (table above).

### 2. Extended introduction / background
- Full prior-art treatment and design rationale, condensed from `docs/EXTENDED_INTRODUCTION.md`. (Supplementary if the venue has no such section.)

### 3. Methods
- GWAS inputs: FinnGen R11 `G6_ALZHEIMER` (11,755 cases / 441,978 controls), `G6_PARKINSON` (5,150 / 448,583); deterministic QC.
- Cell-type programs: MapMyCells top-100 ranked markers per taxonomy node, protein-coding only; testability filters (30/31 superclusters, 436/461 clusters pass ≥10-gene filter).
- Enrichment test: SNP→gene ±10 kb, mean SNP chi-square gene statistic, per-program OLS with log-SNP-count covariate; BH-FDR α = 0.05 per family.
- Robustness: top-100 vs top-200 marker sensitivity.
- Validation: synthetic GWAS with planted enrichment, recovered by both S-LDSC and MAGMA paths (`tests/`).
- Stated as pending at submission: official S-LDSC run with external LD reference (scripted, `scripts/run_ldsc.sh`), and replication on EBI-staged GWAS — framed as pre-registered next steps, not hidden gaps.

### 4. Results

Every table below already exists in captioned form in `docs/RESULTS.md`; Figures 1–3 are committed as reviewable SVGs and regenerate deterministically with the scripts tracked in issue #7. Figure 4 remains blocked by the uncommitted PD cluster table (issue #10).

| # | Content | Source |
|---|---|---|
| Table 1 | Run overview: traits, GWAS sources, QC counts, genes/programs tested, FDR hit counts per level | `docs/RESULTS.md` Table 1 |
| Table 2 | AD supercluster hits (4/30): OPC, Vascular, Microglia, Astrocyte | `docs/RESULTS.md` Table 2 |
| Table 3 | AD cluster hits (9/436), glia-dominated; top `Mgl_9` (p = 2.2e-21) | `docs/RESULTS.md` Table 3 |
| Table 4 | PD summary: 0 class / 0 supercluster / 11 cluster hits; top `CA4_197` (p = 5.3e-13) | `docs/RESULTS.md` Table 4 |
| Table 5 | Resolution gain: fine vs coarse hit counts and minimum p per trait — **the headline table** | `docs/RESULTS.md` Table 5 |
| Table 6 | Robustness: Spearman ρ = 0.81 (AD) / 0.86 (PD), top-100 vs top-200 z-scores | `docs/RESULTS.md` Table 6 |
| Fig. 1 | Enrichment atlas overview across trait × level (hit count + minimum p heatmap) | `figures/figure1_atlas_overview.svg` |
| Fig. 2 | Resolution-gain paired plot (supercluster vs cluster, per trait) | `figures/figure2_resolution_gain.svg` |
| Fig. 3 | AD cluster-level forest plot with 95% intervals | `figures/figure3_AD_FinnGen_R11_cluster_forest.svg` |
| Fig. 4 | PD cluster-level forest plot — **blocked until the PD cluster table is committed (issue #10)** | `figures/figure4_PD_FinnGen_R11_cluster_forest.svg` (pending) |

### 5. Discussion
- Methodological claim: fine taxonomy resolution is a requirement, not a refinement (PD signal exists *only* at cluster level).
- AD glial dominance (four microglial clusters + OPC + pericyte + astrocytes) is concordant with the literature — a positive control on the method.
- PD cluster hits as hypothesis-generating; no mechanistic claims (association only).
- Limitations stated plainly: in-repo MAGMA-style pipeline pending official S-LDSC concordance; single GWAS source per trait; locus-to-gene follow-up deferred to the series' later papers.

### 6. Data and code availability
- All inputs md5-pinned; run metadata in `reports/RUN_METADATA.json`; frozen-manifest re-run procedure and its verification log (issue #8).

### 7. References (core)
- Siletti et al. 2023 (ABC Atlas taxonomy); S-LDSC, MAGMA, sc-linker method papers; FinnGen R11 release.

---

## Pre-submission checklist (tracked as issues)

- [x] This outline committed (`docs/MANUSCRIPT_OUTLINE.md`) — issue #6
- [ ] Figure scripts producing Fig. 1–4 deterministically — Figures 1–3 committed as SVGs; Figure 4 awaits the PD cluster table (#10) — issue #7
- [ ] PD cluster-level full table committed (removes the known gap in `docs/RESULTS.md` Table 4) — issue #10
- [ ] Clean-room reproducibility re-run log — issue #8
- [ ] Official S-LDSC run with external LD reference — scripted; scheduled
- [ ] Scientific owner review before any submission step
