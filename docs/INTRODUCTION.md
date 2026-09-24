# Introduction — Paper 1: Cell-Type-Resolved Heritability Atlas of Brain Traits

**Series note:** This is **Paper 1 of 4** in the GWAS × brain cell-type series. It is the first paper and does not build on any former paper. Papers 2 (cross-species mapping), 3 (coloc/TWAS prioritization), and 4 (psychiatric extension) all reuse the harmonized GWAS inputs and cell-type gene programs built here.

**Concept figure:** see [figures/concept_figure.md](figures/concept_figure.md) for the canonical Mermaid concept diagram.

*Concept figure: GWAS summary statistics and ABC Atlas cell-type marker programs are combined into gene-level disease scores and tested for enrichment, yielding one heritability-enrichment score per cell type at each level of the cell-type hierarchy. (Mermaid-renderable version: [figures/concept_figure.md](figures/concept_figure.md).)*

## Background

Genome-wide association studies (GWAS) have transformed the genetics of brain disease. For Alzheimer's disease (AD), Parkinson's disease (PD), schizophrenia, depression, and many other neuropsychiatric and neurodegenerative traits, large meta-analyses have mapped dozens to hundreds of genome-wide significant loci [8,12,13]. Yet two related problems remain. First, common-variant GWAS identify associated *loci*, not causal genes, cell types, or mechanisms. Second, most trait-associated variants lie outside protein-coding exons, in noncoding regions whose functional interpretation depends entirely on regulatory annotation [14]. The question this paper asks is the central cell-type question of statistical genetics: *which cell types in the human brain carry the genetic risk for these traits?*

The link between noncoding variants and cell types runs through gene regulation. Consortia such as the Roadmap Epigenomics Consortium have shown that trait-associated variants are systematically enriched in cell-type-specific regulatory elements — enhancers and promoters marked by chromatin state — of the tissues plausibly causal for the trait [14]. Stratified linkage disequilibrium score regression (S-LDSC) formalized this observation into partitioned heritability: GWAS summary statistics are regressed on LD scores computed within genomic annotations, estimating how much each annotation contributes to total SNP heritability [1]. Extending the annotation set to genes specifically expressed in a given tissue or cell type recovered, for many traits, the correct causal tissue — and, in follow-up work, specific neuronal and glial populations for brain traits [2,6,7]. Gene-level methods such as MAGMA provide a complementary route, aggregating SNP associations within gene bodies and testing whether trait-associated genes concentrate in cell-type-specific gene sets [3]. The sc-linker framework goes one step further by connecting single-cell expression programs to disease via epigenomic SNP-to-gene maps, identifying not just cell types but cellular processes carrying disease heritability [4].

Until recently, the limiting input was the cell-type reference itself. Single-cell and single-nucleus RNA and chromatin-accessibility profiling (snRNA-seq, snATAC-seq) now resolve brain cell diversity at unprecedented scale. The Allen Brain Cell (ABC) Atlas whole-human-brain component — approximately three million nuclei sampled across ~100 dissections spanning the forebrain, midbrain, hindbrain, and spinal cord — defines a hierarchical taxonomy of superclusters, clusters, and more than 3,000 cell types [5]. Critically, the taxonomy resolves not only coarse classes (excitatory vs. inhibitory neurons, microglia, astrocytes, oligodendrocytes) but region-specific neuronal subtypes and glial *states*. Each level of this hierarchy defines cell-type-specific gene sets and accessible-chromatin annotations that can be converted directly into S-LDSC annotations and MAGMA gene sets. This makes a whole-hierarchy, fine-resolution heritability atlas computable for the first time.

The logic of this paper is therefore straightforward: if a trait's heritability is disproportionately explained by variants in or near genes specifically expressed (or specifically accessible) in a particular cell type, that cell type is a candidate mediator of genetic risk. Partitioned heritability with S-LDSC quantifies this at the SNP-annotation level [1]; MAGMA tests it at the gene level [3]; sc-linker tests it through enhancer-gene linkage and disease-associated programs [4]. Applying all three to the same cell-type taxonomy and the same harmonized GWAS inputs yields triangulated, method-robust cell-type assignments.

## Prior work and gap

Skene et al. [7] first mapped schizophrenia loci onto specific brain cell types using expression-weighted enrichment, implicating pyramidal cells, medium spiny neurons, and interneurons. Bryois et al. [6] extended this approach across a panel of brain traits, identifying, among other findings, oligodendrocyte signals for multiple sclerosis and neuronal signals for psychiatric traits, using single-cell references of roughly a few hundred types. Finucane et al. established the general S-LDSC framework [1] and the specifically-expressed-gene extension [2]; de Leeuw et al. provided MAGMA [3]; Jagadeesh et al. added sc-linker's epigenomic linkage layer [4].

Three gaps remain. (1) **Resolution:** prior brain-trait enrichment studies used references with tens to a few hundred cell types; the ABC taxonomy offers >3,000, including glial states relevant to AD and PD that coarse labels blur [5]. (2) **Consistency:** prior studies differed in GWAS versions, reference panels, and methods, making results hard to compare across traits; a single harmonized pipeline across methods removes this confound. (3) **Hierarchy:** enrichment at one level (e.g., "microglia") may mask heterogeneity at another (e.g., disease-associated microglial states); no prior work has systematically tested heritability across a full cell-type hierarchy with formal level-aware multiple-testing control.

## Research questions

1. Which human brain cell types, at each level of the ABC Atlas taxonomy, show significant heritability enrichment for AD, PD, and related brain traits?
2. Do AD risk variants concentrate in specific microglial/astrocyte *states* rather than glial types as a whole, and does PD signal extend beyond midbrain dopaminergic neurons?
3. How much additional resolution do S-LDSC, MAGMA, and sc-linker each recover at fine vs. coarse taxonomy levels, and where do they agree?
4. What fraction of trait heritability is attributable to the top enriched cell types?

## Data

| Resource | Content | Scale | Access |
|---|---|---|---|
| Allen Brain Cell Atlas — human whole brain (Siletti et al. 2023) [5] | snRNA-seq taxonomy | ~3M nuclei, >3,000 types | Open (ABC Atlas / CELLxGENE) |
| AD GWAS (Bellenguez et al. 2022) [12] | Summary statistics | ~111k cases / ~677k controls | Open (GWAS Catalog, GCST90027158) |
| PD GWAS (Nalls et al. 2019) [13] | Summary statistics | ~33k cases + proxy cases | Open |
| Schizophrenia GWAS (Trubetskoy et al. 2022) [8] | Summary statistics | ~76k cases / ~244k controls | Open (PGC) |
| Depression GWAS (Wray et al. 2018 [9]; Howard et al. 2019 [10]) | Summary statistics | up to ~246k cases | Open (PGC) |
| Roadmap Epigenomics chromatin maps [14] | Baseline regulatory annotations | 111 reference epigenomes | Open |
| 1000 Genomes EUR reference | LD panel | 489 individuals | Open |

## Methods

GWAS summary statistics are harmonized (build liftover, allele harmonization, munging to LDSC format), with LD Score regression intercepts used to verify that observed inflation reflects polygenicity rather than residual confounding [11]. From the ABC taxonomy we derive, at every hierarchy level, (i) top specifically expressed genes per cell type (for MAGMA gene sets and specifically-expressed-gene S-LDSC annotations [2]) and (ii) cell-type-specific accessible regions from snATAC-seq where available. S-LDSC with the baseline-LD model estimates per-annotation heritability enrichment and coefficient significance [1]; MAGMA gene analysis followed by gene-set regression tests cell-type gene sets conditioning on expression covariates [3]; sc-linker builds enhancer-gene-linked programs and tests their disease association [4]. Multiple testing is controlled with FDR within each trait × level, with a Bonferroni check across the full hierarchy. Method concordance is assessed per cell type, and calibration is validated with positive-control contrasts (e.g., blood traits against immune cell types) and matched negative controls.

## Expected contributions

1. A fine-grained, method-triangulated cell-type heritability atlas for major neurodegenerative and neuropsychiatric traits across the ABC hierarchy.
2. Quantification of resolution gains from atlas-scale cell-type references over legacy coarse labels.
3. Versioned, reusable GWAS inputs and cell-type gene programs that Papers 2–4 consume directly.
4. An open, reproducible S-LDSC/MAGMA/sc-linker pipeline for hierarchical single-cell enrichment.

## Scope and boundary

Enrichment is statistical association, not causal mechanism; functional validation is out of scope. Analyses use European-ancestry GWAS and LD references, limiting cross-ancestry claims. Only open data are used. Colocalization and gene prioritization at individual loci are deferred to Paper 3; cross-species questions to Paper 2.

## References

1. Finucane HK, Bulik-Sullivan B, Gusev A, et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. *Nature Genetics* 47:1228–1235 (2015). doi:10.1038/ng.3404
2. Finucane HK, Reshef YA, Anttila V, et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. *Nature Genetics* 50:621–629 (2018). doi:10.1038/s41588-018-0081-4
3. de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data. *PLoS Computational Biology* 11:e1004219 (2015). doi:10.1371/journal.pcbi.1004219
4. Jagadeesh KA, Dey KK, Montoro DT, et al. Identifying disease-critical cell types and cellular processes by integrating single-cell RNA-sequencing and human genetics. *Nature Genetics* 54:1479–1492 (2022). doi:10.1038/s41588-022-01187-9
5. Siletti K, Hodge R, Mossi Albiach A, et al. Transcriptomic diversity of cell types across the adult human brain. *Science* 382:eadd7046 (2023). doi:10.1126/science.add7046
6. Bryois J, Skene NG, Hansen TF, et al. Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson's disease. *Nature Genetics* 52:482–493 (2020). doi:10.1038/s41588-020-0610-9
7. Skene NG, Bryois J, Bakken TE, et al. Genetic identification of brain cell types underlying schizophrenia. *Nature Genetics* 50:825–833 (2018). doi:10.1038/s41588-018-0129-5
8. Trubetskoy V, Pardiñas AF, Qi T, et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature* 604:502–508 (2022). doi:10.1038/s41586-022-04434-5
9. Wray NR, Ripke S, Mattheisen M, et al. Genome-wide association analyses identify 44 risk variants and refine the genetic architecture of major depression. *Nature Genetics* 50:668–681 (2018). doi:10.1038/s41588-018-0090-3
10. Howard DM, Adams MJ, Clarke TK, et al. Genome-wide meta-analysis of depression identifies 102 independent variants and highlights the importance of the prefrontal brain regions. *Nature Neuroscience* 22:343–352 (2019). doi:10.1038/s41593-018-0326-7
11. Bulik-Sullivan BK, Loh PR, Finucane HK, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics* 47:291–295 (2015). doi:10.1038/ng.3211
12. Bellenguez C, Küçükali F, Jansen IE, et al. New insights into the genetic etiology of Alzheimer's disease and related dementias. *Nature Genetics* 54:412–436 (2022). doi:10.1038/s41588-022-01024-9
13. Nalls MA, Blauwendraat C, Vallerga CL, et al. Identification of novel risk loci, causal insights, and heritable risk for Parkinson's disease: a meta-analysis of genome-wide association studies. *The Lancet Neurology* 18:1091–1102 (2019). doi:10.1016/S1474-4422(19)30320-5
14. Roadmap Epigenomics Consortium. Integrative analysis of 111 reference human epigenomes. *Nature* 518:317–330 (2015). doi:10.1038/nature14248
