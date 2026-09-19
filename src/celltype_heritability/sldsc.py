"""S-LDSC-style stratified heritability enrichment regression.

Implements the core weighted regression of Finucane et al. (2015):

    E[chi2_j] = N * sum_k tau_k * l2(j, k) + 1

with annotation coefficients ``tau_k`` estimated by weighted least squares
and block-jackknife standard errors. Production analyses use the official
``ldsc``; this wrapper exists for synthetic tests and prototyping.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class SLDSCResult:
    """Per-annotation coefficient estimates with jackknife standard errors."""

    coefficients: pd.Series
    standard_errors: pd.Series
    z_scores: pd.Series
    p_values: pd.Series
    n_blocks: int


def _wls(y: np.ndarray, x: np.ndarray, w: np.ndarray) -> np.ndarray:
    sw = np.sqrt(w)
    xw, yw = x * sw[:, None], y * sw
    coef, *_ = np.linalg.lstsq(xw, yw, rcond=None)
    return coef


def run_sldsc_regression(
    chisq: np.ndarray,
    ld_scores: pd.DataFrame,
    n_samples: float,
    n_blocks: int = 50,
) -> SLDSCResult:
    """Fit E[chi2] = 1 + N * sum_k tau_k * l2_k by WLS with jackknife SEs.

    Parameters
    ----------
    chisq : (n_snps,) chi-square association statistics.
    ld_scores : DataFrame (n_snps, n_annots) of partitioned LD scores.
        A free intercept column is added automatically (do not include one).
    n_samples : GWAS sample size.
    n_blocks : number of contiguous blocks for the delete-one jackknife.

    Returns
    -------
    SLDSCResult with per-annotation tau coefficients (per-allele scale).
    """
    from scipy import stats

    chisq = np.asarray(chisq, dtype=float)
    l2 = ld_scores.to_numpy(dtype=float)
    m, k = l2.shape
    if chisq.shape[0] != m:
        raise ValueError("chisq and ld_scores must have the same SNP count")
    x = np.column_stack([np.ones(m), l2]) * np.array([1.0] + [n_samples] * k)
    # S-LDSC-style weights: downweight high-LD / high-chisq SNPs
    w = 1.0 / (1.0 + l2.sum(axis=1)) / np.maximum(chisq, 1.0)

    coef = _wls(chisq, x, w)
    coef[1:] /= n_samples  # report tau on the per-allele scale

    # block jackknife over contiguous SNP blocks
    n_blocks = min(n_blocks, m)
    blocks = np.array_split(np.arange(m), n_blocks)
    jack = np.empty((n_blocks, x.shape[1]))
    for b, idx in enumerate(blocks):
        keep = np.ones(m, dtype=bool)
        keep[idx] = False
        c = _wls(chisq[keep], x[keep], w[keep])
        c[1:] /= n_samples
        jack[b] = c
    se = np.sqrt((n_blocks - 1) / n_blocks * ((jack - jack.mean(axis=0)) ** 2).sum(axis=0))

    names = ["intercept", *ld_scores.columns]
    coef_s = pd.Series(coef, index=names)
    se_s = pd.Series(se, index=names)
    z = coef_s / se_s
    p = pd.Series(2 * stats.norm.sf(np.abs(z)), index=names)
    return SLDSCResult(coef_s, se_s, z, p, n_blocks)


def enrichment(result: SLDSCResult, prop_annotated: pd.Series) -> pd.Series:
    """Convert taus to per-SNP heritability enrichment: tau / p_k + 1.

    ``prop_annotated`` gives the genome-wide fraction of SNPs annotated to
    each category (index matching annotation names).
    """
    tau = result.coefficients.drop("intercept")
    m_frac = prop_annotated.reindex(tau.index).to_numpy(dtype=float)
    return 1.0 + tau.to_numpy() / m_frac
