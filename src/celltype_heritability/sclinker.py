"""sc-linker-style disease-relevance scoring of cell-type programs (skeleton).

Follows Jagadeesh et al. (2022): cell-type expression programs are combined
with gene-level disease association to score each cell type's relevance.
This is a minimal, deterministic skeleton for synthetic testing; production
analyses use the official sc-linker pipeline.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def score_programs(
    expression: pd.DataFrame,
    program_genes: dict[str, list[str]],
) -> pd.DataFrame:
    """Mean scaled expression of each gene program per cell type.

    Parameters
    ----------
    expression : genes x cell-types DataFrame of (normalized) expression.
    program_genes : mapping of program name -> member genes.

    Returns
    -------
    DataFrame (cell types x programs) of mean z-scored expression.
    """
    z = expression.sub(expression.mean(axis=1), axis=0).div(
        expression.std(axis=1).replace(0, np.nan), axis=0
    )
    return pd.DataFrame(
        {
            name: z.loc[z.index.intersection(members)].mean(axis=0)
            for name, members in program_genes.items()
        },
        index=expression.columns,
    )


def disease_relevance_scores(
    program_scores: pd.DataFrame,
    gene_disease_scores: pd.Series,
    program_genes: dict[str, list[str]],
) -> pd.Series:
    """Weight each program by the disease association of its member genes.

    Combines the per-cell-type program scores with gene-level disease
    association statistics (e.g. MAGMA gene z-scores) into one relevance
    score per cell type. Deterministic given its inputs.
    """
    weights = {
        name: gene_disease_scores.reindex(
            gene_disease_scores.index.intersection(members)
        ).mean()
        for name, members in program_genes.items()
    }
    w = pd.Series(weights).fillna(0.0)
    return program_scores.mul(w, axis=1).sum(axis=1)
