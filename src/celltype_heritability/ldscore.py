"""Partitioned LD-score computation.

.. note::
   Production analyses in this project use the official ``ldsc`` software
   (Bulik-Sullivan et al. 2015) against the 1000 Genomes reference panel.
   This module provides a lightweight, dependency-free toy implementation
   used **only** for synthetic tests of the statistics plumbing; it is not a
   replacement for ``ldsc`` on real data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_ld_scores(
    genotypes: np.ndarray,
    snp_positions: np.ndarray,
    annotations: pd.DataFrame | None = None,
    window_bp: int = 1_000_000,
) -> pd.DataFrame:
    """Toy in-window partitioned LD scores.

    For each SNP j and annotation column k::

        l2(j, k) = sum_i r_ij^2 * a_k(i)   over SNPs i within ``window_bp`` of j

    Parameters
    ----------
    genotypes : (n_samples, n_snps) array of genotype dosages.
    snp_positions : (n_snps,) base-pair positions (single chromosome toy model).
    annotations : optional DataFrame (n_snps, n_annots) of 0/1 annotations.
        If None, a single "base" annotation of all ones is used.
    window_bp : flank, in bp, defining the LD window around each SNP.

    Returns
    -------
    DataFrame (n_snps, n_annots) of partitioned LD scores.
    """
    g = np.asarray(genotypes, dtype=float)
    n, m = g.shape
    pos = np.asarray(snp_positions, dtype=float)
    if pos.shape[0] != m:
        raise ValueError("snp_positions length must match genotype columns")
    g = (g - g.mean(axis=0)) / (g.std(axis=0) + 1e-12)
    # pairwise r^2 in one go (fine for toy-sized problems)
    r2 = (g.T @ g / n) ** 2
    in_window = np.abs(pos[:, None] - pos[None, :]) <= window_bp
    r2 = r2 * in_window

    if annotations is None:
        annotations = pd.DataFrame({"base": np.ones(m)}, index=np.arange(m))
    a = annotations.to_numpy(dtype=float)
    l2 = r2 @ a
    return pd.DataFrame(l2, index=annotations.index, columns=annotations.columns)
