"""Build binary cell-type SNP annotations from cell-type marker/peak tables.

The Allen Brain Cell Atlas provides cell-type marker genes and (for ATAC
modalities) cell-type-specific peaks. Both reduce to genomic intervals that
can be turned into binary SNP annotations ("is this SNP in/near a cell-type
marker interval?") for partitioned heritability analyses.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

BED_COLUMNS = ["chrom", "start", "end", "name", "score"]


def table_to_bed(
    table: pd.DataFrame,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    name_col: str | None = None,
    score_col: str | None = None,
) -> pd.DataFrame:
    """Normalize a marker/peak table to a BED-like DataFrame.

    Parameters
    ----------
    table : DataFrame with at least chromosome/start/end columns.
    chrom_col, start_col, end_col : column names carrying the interval.
    name_col : optional column with interval names (e.g. gene symbols).
    score_col : optional column with scores (e.g. marker specificity).

    Returns
    -------
    DataFrame with columns ``BED_COLUMNS`` sorted by (chrom, start).
    """
    missing = {chrom_col, start_col, end_col} - set(table.columns)
    if missing:
        raise ValueError(f"table is missing required columns: {sorted(missing)}")
    bed = pd.DataFrame(
        {
            "chrom": table[chrom_col].astype(str),
            "start": table[start_col].astype(int),
            "end": table[end_col].astype(int),
            "name": table[name_col] if name_col else ".",
            "score": table[score_col] if score_col else 0.0,
        }
    )
    if (bed["end"] <= bed["start"]).any():
        raise ValueError("all intervals must satisfy start < end")
    return bed.sort_values(["chrom", "start"]).reset_index(drop=True)


def add_flank(bed: pd.DataFrame, window: int) -> pd.DataFrame:
    """Extend each interval by ``window`` bp on both sides (clipped at 0)."""
    if window < 0:
        raise ValueError("window must be non-negative")
    out = bed.copy()
    out["start"] = (out["start"] - window).clip(lower=0)
    out["end"] = out["end"] + window
    return out


def binarize_annotation(bed: pd.DataFrame, snps: pd.DataFrame) -> np.ndarray:
    """Return a 0/1 vector: 1 if SNP overlaps any interval in ``bed``.

    ``snps`` must have columns ``chrom`` and ``pos``.
    """
    if not {"chrom", "pos"} <= set(snps.columns):
        raise ValueError("snps must have 'chrom' and 'pos' columns")
    flags = np.zeros(len(snps), dtype=int)
    for chrom, idx in snps.groupby("chrom", sort=False).groups.items():
        chrom_bed = bed[bed["chrom"] == str(chrom)]
        if chrom_bed.empty:
            continue
        starts = chrom_bed["start"].to_numpy()
        ends = chrom_bed["end"].to_numpy()
        pos = snps.loc[idx, "pos"].to_numpy()
        # SNP overlaps if any interval has start <= pos < end
        hit = np.zeros(len(pos), dtype=bool)
        order = np.argsort(starts)
        starts, ends = starts[order], ends[order]
        for s, e in zip(starts, ends):
            hit |= (pos >= s) & (pos < e)
        flags[snps.index.get_indexer(idx)] = hit.astype(int)
    return flags


def build_celltype_annotations(
    beds: dict[str, pd.DataFrame],
    snps: pd.DataFrame,
    window: int = 0,
) -> pd.DataFrame:
    """Build one binary annotation column per cell type.

    Parameters
    ----------
    beds : mapping of cell-type name -> BED-like DataFrame.
    snps : DataFrame with ``chrom`` and ``pos`` columns (SNP universe).
    window : flanking window added to every interval before overlap.

    Returns
    -------
    DataFrame indexed like ``snps`` with one 0/1 column per cell type.
    """
    return pd.DataFrame(
        {
            ct: binarize_annotation(add_flank(bed, window), snps)
            for ct, bed in beds.items()
        },
        index=snps.index,
    )
