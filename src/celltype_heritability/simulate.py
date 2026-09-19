"""Synthetic GWAS summary statistics and annotations with planted enrichment.

Generates a toy genome where SNPs annotated to one designated cell type
carry excess per-SNP heritability, so that both the S-LDSC and MAGMA paths
must recover the planted category.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class SyntheticDataset:
    snps: pd.DataFrame  # chrom, pos
    genotypes: np.ndarray  # n_ref x n_snps toy reference panel
    annotations: pd.DataFrame  # n_snps x n_celltypes binary
    z: np.ndarray  # GWAS z-scores with planted enrichment
    genes: pd.DataFrame  # gene, chrom, start, end
    gene_programs: dict[str, list[str]]  # cell type -> marker genes
    enriched_celltype: str


def simulate_gwas(
    n_snps: int = 4000,
    n_ref: int = 200,
    n_celltypes: int = 4,
    n_genes: int = 400,
    markers_per_type: int = 25,
    sample_size: float = 50_000.0,
    h2_total: float = 0.3,
    enriched_fraction: float = 0.5,
    enriched_celltype: str | None = None,
    seed: int = 0,
) -> SyntheticDataset:
    """Simulate summary statistics with planted cell-type heritability.

    ``enriched_fraction`` of total heritability is assigned uniformly to SNPs
    annotated to ``enriched_celltype`` (default: type 0); the remainder is
    spread uniformly across all SNPs (baseline), so the enriched category
    shows elevated per-SNP heritability relative to the others.
    """
    rng = np.random.default_rng(seed)
    ct_names = [f"celltype_{i}" for i in range(n_celltypes)]
    if enriched_celltype is None:
        enriched_celltype = ct_names[0]

    # toy genome layout: one chromosome, SNPs every 1kb, genes every 10kb
    pos = np.arange(n_snps) * 1_000
    snps = pd.DataFrame({"chrom": "1", "pos": pos})
    genes = pd.DataFrame(
        {
            "gene": [f"gene_{i}" for i in range(n_genes)],
            "chrom": "1",
            "start": np.arange(n_genes) * (n_snps // n_genes) * 1_000,
        }
    )
    genes["end"] = genes["start"] + (n_snps // n_genes) * 1_000

    # sparse genotype panel (toy LD: AR(1)-like correlation along the genome)
    base = rng.normal(size=(n_ref, n_snps))
    genotypes = base.copy()
    for j in range(1, n_snps):
        genotypes[:, j] = 0.3 * genotypes[:, j - 1] + np.sqrt(1 - 0.09) * base[:, j]

    # cell-type programs: each type owns a block of marker genes
    gene_programs = {
        ct: list(genes["gene"][i * markers_per_type : (i + 1) * markers_per_type])
        for i, ct in enumerate(ct_names)
    }

    # binary SNP annotations: SNP is annotated if inside a marker gene
    from .annotations import build_celltype_annotations, table_to_bed

    beds = {
        ct: table_to_bed(genes[genes["gene"].isin(members)], name_col="gene")
        for ct, members in gene_programs.items()
    }
    annot = build_celltype_annotations(beds, snps, window=2_000)

    # per-SNP heritability with planted enrichment
    h2 = np.full(n_snps, (1 - enriched_fraction) * h2_total / n_snps)
    mask = annot[enriched_celltype].to_numpy().astype(bool)
    h2[mask] += enriched_fraction * h2_total / mask.sum()

    # z-scores: E[z^2] = 1 + N * h2_j  (no LD-induced inflation in the toy)
    z = rng.normal(size=n_snps) * np.sqrt(1.0 + sample_size * h2)
    return SyntheticDataset(snps, genotypes, annot, z, genes, gene_programs, enriched_celltype)
