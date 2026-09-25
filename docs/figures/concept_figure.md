# Concept Figure — Cell-Type-Resolved Heritability Atlas

> **Note.** The rendered concept figure is [`01-concept-schematic.png`](01-concept-schematic.png) in this folder — a generated scientific illustration. This file keeps the Mermaid source of the same diagram so it remains editable and re-renderable.

**Caption:** Concept figure — the pipeline: harmonized GWAS summary statistics and Allen Brain Cell Atlas cell-type marker programs are combined into gene-level disease scores and tested for enrichment, yielding one score per cell type at every level of the cell-type hierarchy.

```mermaid
flowchart LR
    subgraph Inputs
        A["GWAS summary statistics<br/>Alzheimer's + Parkinson's<br/>(FinnGen R11, ~21M variants)"]
        B["Allen Brain Cell Atlas<br/>~3M nuclei, >3,000 cell types<br/>class / supercluster / cluster"]
    end
    subgraph Processing
        C["Harmonized gene-level<br/>disease scores<br/>(18,551 protein-coding genes)"]
        D["Cell-type marker gene programs<br/>(top-100 markers per cell type)"]
    end
    E["Enrichment regression<br/>one z-score per cell type<br/>+ FDR control per trait x level"]
    subgraph Outputs
        F["Alzheimer's disease:<br/>OPC z=5.73, Vascular z=5.51,<br/>Microglia z=3.84, Astrocyte z=3.52"]
        G["Parkinson's disease:<br/>11 fine clusters FDR<0.05<br/>(resolution gain 0 -> 11)"]
    end
    A --> C
    B --> D
    C --> E
    D --> E
    E --> F
    E --> G
```
