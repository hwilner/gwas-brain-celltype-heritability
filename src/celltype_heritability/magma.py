"""MAGMA-style gene-level association testing from GWAS summary statistics.

Simplified implementation of de Leeuw et al. (2015): map SNPs to genes with
a flanking window, aggregate SNP chi-squares into gene test statistics, and
test cell-type programs by regressing gene statistics on program membership.
Synthetic-data development only; production runs use official MAGMA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


def map_snps_to_genes(
    snps: pd.DataFrame,
    genes: pd.DataFrame,
    window: int = 10_000,
) -> dict[str, list[int]]:
    """Map SNPs to genes whose (flanked) span contains the SNP position.

    Parameters
    ----------
    snps : DataFrame with columns ``chrom`` and ``pos``.
    genes : DataFrame with columns ``gene``, ``chrom``, ``start``, ``end``.
    window : flank added to both gene boundaries (bp).

    Returns
    -------
    Mapping of gene name -> list of SNP row indices.

    Notes
    -----
    Uses a per-chromosome sorted-position sweep (binary search per gene),
    so it scales to genome-wide SNP sets (millions of SNPs x ~20k genes).
    """
    if not {"chrom", "pos"} <= set(snps.columns):
        raise ValueError("snps must have 'chrom' and 'pos' columns")
    for col in ("gene", "chrom", "start", "end"):
        if col not in genes.columns:
            raise ValueError(f"genes is missing required column {col!r}")

    out: dict[str, list[int]] = {}
    chrom_pos: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, idx in snps.groupby("chrom", sort=False).groups.items():
        idx = np.asarray(idx)
        order = np.argsort(snps["pos"].to_numpy()[snps.index.get_indexer(idx)])
        chrom_pos[str(chrom)] = (
            idx[order],
            snps["pos"].to_numpy()[snps.index.get_indexer(idx)][order],
        )
    for gene, grp in genes.groupby("gene", sort=False):
        hits: list[int] = []
        for _, row in grp.iterrows():
            key = str(row["chrom"])
            if key not in chrom_pos:
                continue
            idx, pos = chrom_pos[key]
            lo = np.searchsorted(pos, row["start"] - window, side="left")
            hi = np.searchsorted(pos, row["end"] + window, side="right")
            hits.extend(idx[lo:hi].tolist())
        out[str(gene)] = sorted(set(hits))
    return out


def gene_test_statistics(z: np.ndarray, snp_map: dict[str, list[int]]) -> pd.Series:
    """Mean SNP chi-square per gene (MAGMA's SNP-wise mean model, no LD)."""
    z = np.asarray(z, dtype=float)
    z2 = z**2
    return pd.Series(
        {g: z2[idx].mean() for g, idx in snp_map.items() if len(idx) > 0}
    )


@dataclass
class MagmaEnrichmentResult:
    coefficient: float
    standard_error: float
    z_score: float
    p_value: float
    n_genes: int


def program_enrichment(
    gene_stats: pd.Series,
    gene_sets: dict[str, set[str]] | dict[str, list[str]],
    gene_sizes: pd.Series | None = None,
) -> pd.DataFrame:
    """Regress gene test statistics on cell-type program membership.

    Model: ``gene_stat ~ 1 + log(n_snps) + 1[gene in program]`` fit by OLS.
    The program coefficient > 0 indicates enrichment of association signal in
    the cell type's marker genes.
    """
    from scipy import stats

    genes = gene_stats.index
    if gene_sizes is None:
        gene_sizes = pd.Series(1.0, index=genes)
    rows = []
    for ct, members in gene_sets.items():
        in_prog = genes.isin(list(members)).astype(float)
        x = np.column_stack(
            [np.ones(len(genes)), np.log(gene_sizes.reindex(genes).fillna(1.0)), in_prog]
        )
        y = gene_stats.to_numpy(dtype=float)
        coef, *_ = np.linalg.lstsq(x, y, rcond=None)
        resid = y - x @ coef
        dof = max(len(y) - x.shape[1], 1)
        sigma2 = resid @ resid / dof
        cov = sigma2 * np.linalg.inv(x.T @ x)
        se = np.sqrt(cov[2, 2])
        z = coef[2] / se
        rows.append(
            MagmaEnrichmentResult(
                coefficient=coef[2],
                standard_error=se,
                z_score=z,
                p_value=2 * stats.norm.sf(abs(z)),
                n_genes=len(genes),
            )
        )
    return pd.DataFrame(rows, index=list(gene_sets))


def bh_fdr(p_values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR-adjusted q-values (deterministic)."""
    p = p_values.to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return pd.Series(out, index=p_values.index)
