"""Enrichment atlas figure scripts for Paper 1 (issue #7).

Reads only committed, versioned inputs from ``reports/`` and regenerates
every results figure in ``figures/`` deterministically:

- one volcano plot per committed (trait x level) enrichment table
  (``magma_enrichment_<TRAIT>_<LEVEL>.csv``),
- one top-programs bar chart per committed (trait x level) table,
- the resolution-gain figure from ``resolution_gain.csv``.

The PD cluster-level table is not yet committed (issue #10), so its panels
are skipped with a logged warning until that gap is closed. Figures only;
no numbers are recomputed here beyond plotting transforms of the committed
tables (e.g. ``-log10(p)``), and nothing is simulated or estimated.

Usage::

    python -m celltype_heritability.figures --reports-dir reports --out-dir figures
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless, deterministic backend; must precede pyplot

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

FDR_ALPHA = 0.05
"""Pre-registered Benjamini-Hochberg threshold per (trait x level) family."""

LEVELS = ("class", "supercluster", "cluster")
"""Hierarchy levels in the pre-registered testing grid."""

REQUIRED_COLUMNS = (
    "program",
    "coefficient",
    "standard_error",
    "p_value",
    "q_value",
)
"""Columns every committed enrichment table must provide."""

# Fixed figure style so re-runs produce byte-identical PNGs (metadata=None
# in the savefig calls removes the version-stamped PNG text chunk).
plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def enrichment_table_path(reports_dir: Path, trait: str, level: str) -> Path:
    """Path of the committed enrichment table for one (trait, level) cell.

    Args:
        reports_dir: Directory holding the committed ``magma_enrichment_*``
            tables (conventionally ``reports/``).
        trait: Trait tag as used in filenames, e.g. ``AD_FinnGen_R11``.
        level: Hierarchy level, one of :data:`LEVELS`.

    Returns:
        Path to ``magma_enrichment_<trait>_<level>.csv``.
    """
    return Path(reports_dir) / f"magma_enrichment_{trait}_{level}.csv"


def load_enrichment_table(path: Path) -> pd.DataFrame:
    """Load one committed enrichment table and validate its schema.

    Args:
        path: CSV produced by ``scripts/run_magma_enrichment.py``.

    Returns:
        Table sorted by p-value with a ``neg_log10_p`` plotting column.

    Raises:
        ValueError: If a required column is missing.
    """
    df = pd.read_csv(path)
    missing = set(REQUIRED_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    df = df.sort_values("p_value", kind="stable").reset_index(drop=True)
    # Clip at the smallest positive p so p = 0 cannot produce inf.
    floor = np.finfo(float).tiny
    df["neg_log10_p"] = -np.log10(df["p_value"].clip(lower=floor))
    return df


def _save(fig: plt.Figure, out_path: Path) -> Path:
    """Save a figure deterministically (no timestamp/version metadata)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)
    return out_path


def plot_volcano(df: pd.DataFrame, out_path: Path, title: str) -> Path:
    """Volcano plot: enrichment coefficient vs -log10(p), FDR hits in colour.

    Args:
        df: Table from :func:`load_enrichment_table`.
        out_path: Destination PNG.
        title: Axes title, naming trait and level.

    Returns:
        The written path.
    """
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    hit = df["q_value"] < FDR_ALPHA
    ax.scatter(
        df.loc[~hit, "coefficient"],
        df.loc[~hit, "neg_log10_p"],
        s=8,
        c="0.6",
        linewidths=0,
        label=f"q >= {FDR_ALPHA}",
    )
    ax.scatter(
        df.loc[hit, "coefficient"],
        df.loc[hit, "neg_log10_p"],
        s=14,
        c="crimson",
        linewidths=0,
        label=f"q < {FDR_ALPHA}",
    )
    for _, row in df.loc[hit].head(5).iterrows():
        ax.annotate(
            str(row["program"]),
            (row["coefficient"], row["neg_log10_p"]),
            fontsize=6,
            xytext=(3, 3),
            textcoords="offset points",
        )
    ax.axvline(0.0, color="0.3", linewidth=0.6)
    ax.margins(x=0.12)  # keep top-hit annotations inside the axes
    ax.set_xlabel("Program coefficient (beta)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    return _save(fig, out_path)


def plot_top_programs(df: pd.DataFrame, out_path: Path, title: str, n: int = 10) -> Path:
    """Horizontal bar chart of the n strongest programs (beta with SE bars).

    FDR-significant programs are coloured; the rest are grey. Args and
    return are as in :func:`plot_volcano`; ``n`` caps the panel height.
    """
    top = df.head(n).iloc[::-1]  # strongest on top of a horizontal bar chart
    hit = top["q_value"] < FDR_ALPHA
    colors = np.where(hit, "crimson", "0.6")
    fig, ax = plt.subplots(figsize=(5.6, 0.32 * len(top) + 1.2))
    ax.barh(
        top["program"].astype(str),
        top["coefficient"],
        xerr=top["standard_error"],
        color=colors,
        height=0.7,
        error_kw={"linewidth": 0.8},
    )
    ax.axvline(0.0, color="0.3", linewidth=0.6)
    ax.set_xlabel("Program coefficient (beta +/- SE)")
    ax.set_title(title)
    return _save(fig, out_path)


def plot_resolution_gain(df: pd.DataFrame, out_path: Path) -> Path:
    """Paired bars: FDR-significant programs at supercluster vs cluster level.

    Args:
        df: ``reports/resolution_gain.csv`` with columns ``trait``,
            ``n_sig_supercluster``, ``n_sig_cluster``.
        out_path: Destination PNG.

    Returns:
        The written path.
    """
    traits = df["trait"].astype(str).tolist()
    x = np.arange(len(traits))
    width = 0.36
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.bar(x - width / 2, df["n_sig_supercluster"], width, label="supercluster", color="0.6")
    ax.bar(x + width / 2, df["n_sig_cluster"], width, label="cluster", color="crimson")
    for xi, (_, row) in zip(x, df.iterrows()):
        ax.text(xi - width / 2, row["n_sig_supercluster"] + 0.15, str(int(row["n_sig_supercluster"])), ha="center", fontsize=8)
        ax.text(xi + width / 2, row["n_sig_cluster"] + 0.15, str(int(row["n_sig_cluster"])), ha="center", fontsize=8)
    ax.set_xticks(x, traits)
    ax.set_ylabel(f"Programs with FDR q < {FDR_ALPHA}")
    ax.set_title("Resolution gain: fine vs coarse taxonomy")
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, out_path)


def discover_traits(reports_dir: Path) -> list[str]:
    """Trait tags of committed enrichment tables, in stable sorted order."""
    tags = set()
    for path in Path(reports_dir).glob("magma_enrichment_*.csv"):
        m = re.fullmatch(r"magma_enrichment_(.+)_(class|supercluster|cluster)\.csv", path.name)
        if m:
            tags.add(m.group(1))
    return sorted(tags)


def build_atlas_figures(reports_dir: Path, out_dir: Path) -> list[Path]:
    """Regenerate every enrichment atlas figure from committed tables.

    Missing (trait x level) tables are skipped with a warning rather than
    treated as an error, so the known PD cluster-level gap (issue #10) does
    not block the committed cells.

    Args:
        reports_dir: Directory of committed ``reports/`` tables.
        out_dir: Figure destination (conventionally ``figures/``).

    Returns:
        Sorted list of written PNG paths.

    Raises:
        FileNotFoundError: If ``reports_dir`` or ``resolution_gain.csv``
            does not exist.
    """
    reports_dir = Path(reports_dir)
    out_dir = Path(out_dir)
    if not reports_dir.is_dir():
        raise FileNotFoundError(f"reports directory not found: {reports_dir}")
    gain_path = reports_dir / "resolution_gain.csv"
    if not gain_path.is_file():
        raise FileNotFoundError(f"resolution-gain table not found: {gain_path}")

    written: list[Path] = []
    for trait in discover_traits(reports_dir):
        for level in LEVELS:
            path = enrichment_table_path(reports_dir, trait, level)
            if not path.is_file():
                logger.warning("skipping %s x %s: %s not committed", trait, level, path.name)
                continue
            df = load_enrichment_table(path)
            title = f"{trait} - {level} level"
            written.append(plot_volcano(df, out_dir / f"volcano_{trait}_{level}.png", title))
            written.append(plot_top_programs(df, out_dir / f"top_programs_{trait}_{level}.png", title))
    gain = pd.read_csv(gain_path)
    written.append(plot_resolution_gain(gain, out_dir / "resolution_gain.png"))
    return sorted(written)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: regenerate all atlas figures from ``reports/``."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reports-dir", type=Path, default=Path("reports"))
    parser.add_argument("--out-dir", type=Path, default=Path("figures"))
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    written = build_atlas_figures(args.reports_dir, args.out_dir)
    for path in written:
        logger.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
