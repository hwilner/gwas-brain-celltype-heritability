"""Enrichment atlas figure scripts for Paper 1 (issue #7).

Reads only committed, versioned inputs from ``reports/`` and regenerates the
Paper 1 figure package deterministically:

- Figure 1, an atlas overview heatmap built from ``resolution_gain.csv``;
- Figure 2, the paired supercluster-vs-cluster resolution-gain plot;
- Figure 3 onward, forest plots (beta with 95% confidence intervals) for each
  committed cluster-level enrichment table;
- diagnostic volcano and top-program panels for every committed
  (trait x level) enrichment table.

Every figure is written as both PNG (for local inspection) and SVG (the
committed, reviewable text format). The PD cluster-level table is not yet
committed (issue #10), so its panels are skipped with a logged warning until
that gap is closed. Figures only; no numbers are recomputed here beyond
plotting transforms of the committed tables (e.g. ``-log10(p)``), and nothing
is simulated or estimated.

Usage::

    python -m celltype_heritability.figures --reports-dir reports --out-dir figures
"""

from __future__ import annotations

import argparse
import logging
import re
from html import escape
from pathlib import Path
from typing import Callable

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

FIGURE_FORMATS = ("png", "svg")
"""Output formats written by default: local PNG plus committed SVG."""

REQUIRED_COLUMNS = (
    "program",
    "coefficient",
    "standard_error",
    "p_value",
    "q_value",
)
"""Columns every committed enrichment table must provide."""

REQUIRED_GAIN_COLUMNS = (
    "trait",
    "n_sig_supercluster",
    "n_sig_cluster",
    "min_p_supercluster",
    "min_p_cluster",
)
"""Columns every resolution-gain table must provide."""

# Fixed figure style so re-runs produce byte-identical outputs in a fixed
# environment. ``svg.fonttype=none`` keeps committed SVG text selectable.
plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        # Matplotlib otherwise salts SVG element IDs with a per-process random
        # value, which makes byte-identical regeneration impossible.
        "svg.hashsalt": "celltype-heritability",
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


def load_resolution_gain(path: Path) -> pd.DataFrame:
    """Load and validate ``reports/resolution_gain.csv``.

    Args:
        path: Resolution-gain CSV produced by
            ``scripts/run_magma_enrichment.py``.

    Returns:
        The validated table in committed row order.

    Raises:
        ValueError: If a required column is missing or a minimum p-value is
            outside ``(0, 1]``.
    """
    df = pd.read_csv(path)
    missing = set(REQUIRED_GAIN_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    for col in ("min_p_supercluster", "min_p_cluster"):
        values = df[col].to_numpy(dtype=float)
        if not np.all(np.isfinite(values)) or not np.all((values > 0) & (values <= 1)):
            raise ValueError(f"{path} has invalid minimum p-values in {col!r}")
    return df


def _save(fig: plt.Figure, out_path: Path, title: str | None = None) -> Path:
    """Save a figure deterministically (no timestamp metadata).

    PNG accepts a ``Software`` metadata key; SVG uses Dublin Core keys. Setting
    SVG ``Date`` to ``None`` removes the wall-clock timestamp, and pinning
    ``Creator`` avoids embedding a matplotlib version string in committed
    figures.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.suffix == ".png":
        metadata = {"Software": None}
    elif out_path.suffix == ".svg":
        metadata = {
            "Creator": "celltype-heritability figures.py",
            "Date": None,
            "Title": title or out_path.stem,
        }
    else:
        raise ValueError(f"unsupported figure format: {out_path.suffix}")
    fig.savefig(out_path, bbox_inches="tight", metadata=metadata)
    plt.close(fig)
    return out_path


def _write_svg(
    out_path: Path,
    width: int,
    height: int,
    title: str,
    body: list[str],
) -> Path:
    """Write a small, deterministic vector SVG for a publication figure."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    content = "\n".join(
        [
            '<?xml version="1.0" encoding="UTF-8"?>',
            (
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
                f'height="{height}" viewBox="0 0 {width} {height}" role="img" '
                f'aria-label="{escape(title)}">'
            ),
            f"<title>{escape(title)}</title>",
            '<rect width="100%" height="100%" fill="white"/>',
            *body,
            "</svg>",
            "",
        ]
    )
    out_path.write_text(content, encoding="utf-8")
    return out_path


def _svg_text(x: float, y: float, text: object, size: int = 11, anchor: str = "start", fill: str = "#222") -> str:
    """One SVG text element with the figure package's fixed font stack."""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" text-anchor="{anchor}" fill="{fill}">{escape(str(text))}</text>'
    )


def _heat_color(value: float, vmin: float, vmax: float) -> str:
    """Map a scalar to a deterministic light-to-dark blue hex colour."""
    t = 0.0 if vmax == vmin else (value - vmin) / (vmax - vmin)
    lo = np.array([247, 251, 255])
    hi = np.array([8, 48, 107])
    rgb = lo + t * (hi - lo)
    return "#" + "".join(f"{int(round(c)):02x}" for c in rgb)


def plot_atlas_overview_svg(df: pd.DataFrame, out_path: Path) -> Path:
    """Write Figure 1 as a compact vector SVG from resolution-gain values."""
    levels = ("supercluster", "cluster")
    count_cols = {"supercluster": "n_sig_supercluster", "cluster": "n_sig_cluster"}
    p_cols = {"supercluster": "min_p_supercluster", "cluster": "min_p_cluster"}
    scores = -np.log10(df[[p_cols[level] for level in levels]].to_numpy(dtype=float))
    vmin, vmax = float(scores.min()), float(scores.max())
    x0, y0, cell_w, cell_h = 250.0, 65.0, 150.0, 62.0
    body = [_svg_text(320, 28, "Enrichment atlas overview", 16, "middle")]
    for j, level in enumerate(levels):
        body.append(_svg_text(x0 + j * cell_w + cell_w / 2, y0 - 12, level, 11, "middle"))
    for i, (_, row) in enumerate(df.iterrows()):
        y = y0 + i * cell_h
        body.append(_svg_text(x0 - 12, y + cell_h / 2 + 4, row["trait"], 11, "end"))
        for j, level in enumerate(levels):
            x = x0 + j * cell_w
            score = float(scores[i, j])
            fill = _heat_color(score, vmin, vmax)
            t = 0.0 if vmax == vmin else (score - vmin) / (vmax - vmin)
            text_fill = "white" if t > 0.55 else "#222"
            body.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{cell_w:.1f}" '
                f'height="{cell_h:.1f}" fill="{fill}" stroke="#333" stroke-width="0.7"/>'
            )
            body.append(
                _svg_text(x + cell_w / 2, y + 25, f"n={int(row[count_cols[level]])}", 11, "middle", text_fill)
            )
            body.append(
                _svg_text(x + cell_w / 2, y + 43, f"p={float(row[p_cols[level]]):.1e}", 10, "middle", text_fill)
            )
    body.append(_svg_text(x0, y0 + len(df) * cell_h + 32, "Colour: -log10(minimum p); labels show FDR hit count and minimum p.", 10))
    return _write_svg(out_path, 640, int(y0 + len(df) * cell_h + 52), "Enrichment atlas overview", body)


def plot_resolution_gain_svg(df: pd.DataFrame, out_path: Path) -> Path:
    """Write Figure 2 as a compact vector SVG from committed hit counts."""
    width, height = 640, 360
    x0, y_base, plot_w, plot_h = 90.0, 285.0, 440.0, 220.0
    max_count = int(max(df["n_sig_supercluster"].max(), df["n_sig_cluster"].max(), 1))
    tick_step = max(1, int(np.ceil(max_count / 5)))
    top = tick_step * int(np.ceil(max_count / tick_step))
    body = [_svg_text(width / 2, 28, "Resolution gain: fine vs coarse taxonomy", 16, "middle")]
    body.append(f'<line x1="{x0}" y1="65" x2="{x0}" y2="{y_base}" stroke="#222"/>')
    body.append(f'<line x1="{x0}" y1="{y_base}" x2="{x0 + plot_w}" y2="{y_base}" stroke="#222"/>')
    for tick in range(0, top + 1, tick_step):
        y = y_base - plot_h * tick / top
        body.append(f'<line x1="{x0 - 5}" y1="{y:.1f}" x2="{x0}" y2="{y:.1f}" stroke="#222"/>')
        body.append(_svg_text(x0 - 10, y + 4, tick, 10, "end"))
    group_w = plot_w / len(df)
    bar_w = min(70.0, group_w / 3)
    for i, (_, row) in enumerate(df.iterrows()):
        center = x0 + group_w * (i + 0.5)
        for offset, col, fill in (
            (-bar_w / 1.8, "n_sig_supercluster", "#777"),
            (bar_w / 1.8, "n_sig_cluster", "#c9184a"),
        ):
            count = int(row[col])
            h = plot_h * count / top
            x = center + offset - bar_w / 2
            y = y_base - h
            body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{fill}"/>')
            body.append(_svg_text(x + bar_w / 2, y - 6, count, 10, "middle"))
        body.append(_svg_text(center, y_base + 22, row["trait"], 10, "middle"))
    body.append(f'<rect x="545" y="75" width="12" height="12" fill="#777"/><text x="562" y="85" font-family="Arial" font-size="10">supercluster</text>')
    body.append(f'<rect x="545" y="96" width="12" height="12" fill="#c9184a"/><text x="562" y="106" font-family="Arial" font-size="10">cluster</text>')
    body.append(_svg_text(25, 175, "FDR q < 0.05 programs", 11, "middle"))
    return _write_svg(out_path, width, height, "Resolution gain: fine vs coarse taxonomy", body)


def plot_cluster_forest_svg(df: pd.DataFrame, out_path: Path, title: str, max_rows: int = 12) -> Path:
    """Write a cluster-level forest plot as a compact vector SVG."""
    hits = df[df["q_value"] < FDR_ALPHA]
    display = (hits if len(hits) else df).head(max_rows).iloc[::-1]
    beta = display["coefficient"].to_numpy(dtype=float)
    lo = beta - 1.96 * display["standard_error"].to_numpy(dtype=float)
    hi = beta + 1.96 * display["standard_error"].to_numpy(dtype=float)
    xmin = min(float(lo.min()), 0.0)
    xmax = max(float(hi.max()), 0.1)
    pad = 0.08 * (xmax - xmin)
    xmin, xmax = xmin - pad, xmax + pad
    width = 780
    height = max(220, 105 + 25 * len(display))
    x0, x1, y0 = 190.0, 625.0, 58.0
    row_h = 25.0
    plot_h = row_h * max(len(display) - 1, 1)

    def sx(value: float) -> float:
        return x0 + (value - xmin) / (xmax - xmin) * (x1 - x0)

    body = [_svg_text(width / 2, 28, title, 16, "middle")]
    body.append(f'<line x1="{x0}" y1="{y0 - 15}" x2="{x0}" y2="{y0 + plot_h + 18}" stroke="#222"/>')
    body.append(f'<line x1="{x0}" y1="{y0 + plot_h + 18}" x2="{x1}" y2="{y0 + plot_h + 18}" stroke="#222"/>')
    for tick in np.linspace(xmin, xmax, 5):
        x = sx(float(tick))
        body.append(f'<line x1="{x:.1f}" y1="{y0 + plot_h + 18}" x2="{x:.1f}" y2="{y0 + plot_h + 23}" stroke="#222"/>')
        body.append(_svg_text(x, y0 + plot_h + 39, f"{tick:.1f}", 10, "middle"))
    zero_x = sx(0.0)
    body.append(f'<line x1="{zero_x:.1f}" y1="{y0 - 15}" x2="{zero_x:.1f}" y2="{y0 + plot_h + 18}" stroke="#666" stroke-width="0.8"/>')
    for i, (_, row) in enumerate(display.iterrows()):
        y = y0 + i * row_h
        colour = "#c9184a" if row["q_value"] < FDR_ALPHA else "#777"
        left, right, centre = sx(float(lo[i])), sx(float(hi[i])), sx(float(beta[i]))
        body.append(f'<line x1="{left:.1f}" y1="{y:.1f}" x2="{right:.1f}" y2="{y:.1f}" stroke="#222"/>')
        body.append(f'<line x1="{left:.1f}" y1="{y - 4:.1f}" x2="{left:.1f}" y2="{y + 4:.1f}" stroke="#222"/>')
        body.append(f'<line x1="{right:.1f}" y1="{y - 4:.1f}" x2="{right:.1f}" y2="{y + 4:.1f}" stroke="#222"/>')
        body.append(f'<circle cx="{centre:.1f}" cy="{y:.1f}" r="4" fill="{colour}"/>')
        body.append(_svg_text(x0 - 12, y + 4, row["program"], 10, "end"))
        body.append(_svg_text(x1 + 12, y + 4, f"q={row['q_value']:.1e}", 10))
    body.append(_svg_text((x0 + x1) / 2, height - 18, "Program coefficient (beta, 95% CI)", 11, "middle"))
    return _write_svg(out_path, width, height, title, body)


def _write_formats(
    plotter: Callable[..., Path],
    *args,
    base_path: Path,
    formats: tuple[str, ...],
    svg_plotter: Callable[..., Path] | None = None,
    **kwargs,
) -> list[Path]:
    """Write one plot in every requested format."""
    written = []
    for fmt in formats:
        out_path = base_path.with_suffix(f".{fmt}")
        writer = svg_plotter if fmt == "svg" and svg_plotter is not None else plotter
        written.append(writer(*args, out_path, **kwargs))
    return written


def plot_volcano(df: pd.DataFrame, out_path: Path, title: str) -> Path:
    """Volcano plot: enrichment coefficient vs -log10(p), FDR hits in colour.

    Args:
        df: Table from :func:`load_enrichment_table`.
        out_path: Destination figure path.
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
    return _save(fig, out_path, title)


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
    return _save(fig, out_path, title)


def plot_atlas_overview(df: pd.DataFrame, out_path: Path) -> Path:
    """Figure 1 heatmap: hit counts and minimum p-values by trait x level.

    The overview deliberately uses ``resolution_gain.csv`` rather than the
    per-program enrichment tables so the PD cluster cell can be shown even
    while its full table remains uncommitted (issue #10).

    Args:
        df: Validated table from :func:`load_resolution_gain`.
        out_path: Destination figure path.

    Returns:
        The written path.
    """
    levels = ("supercluster", "cluster")
    count_cols = {"supercluster": "n_sig_supercluster", "cluster": "n_sig_cluster"}
    p_cols = {"supercluster": "min_p_supercluster", "cluster": "min_p_cluster"}
    p_matrix = df[[p_cols[level] for level in levels]].to_numpy(dtype=float)
    score = -np.log10(np.clip(p_matrix, np.finfo(float).tiny, 1.0))

    title = "Enrichment atlas overview"
    fig, ax = plt.subplots(figsize=(4.8, 1.4 + 0.7 * len(df)))
    image = ax.imshow(score, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(levels)), levels)
    ax.set_yticks(np.arange(len(df)), df["trait"].astype(str))
    ax.set_title(title)
    for i, (_, row) in enumerate(df.iterrows()):
        for j, level in enumerate(levels):
            count = int(row[count_cols[level]])
            p_value = float(row[p_cols[level]])
            ax.text(
                j,
                i,
                f"n={count}\np={p_value:.1e}",
                ha="center",
                va="center",
                color="white" if score[i, j] > score.max() * 0.55 else "black",
                fontsize=8,
            )
    colorbar = fig.colorbar(image, ax=ax, shrink=0.85)
    colorbar.set_label("-log10(minimum p)")
    return _save(fig, out_path, title)


def plot_resolution_gain(df: pd.DataFrame, out_path: Path) -> Path:
    """Paired bars: FDR-significant programs at supercluster vs cluster level.

    Args:
        df: ``reports/resolution_gain.csv`` with columns ``trait``,
            ``n_sig_supercluster``, ``n_sig_cluster``.
        out_path: Destination figure path.

    Returns:
        The written path.
    """
    traits = df["trait"].astype(str).tolist()
    x = np.arange(len(traits))
    width = 0.36
    title = "Resolution gain: fine vs coarse taxonomy"
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.bar(x - width / 2, df["n_sig_supercluster"], width, label="supercluster", color="0.6")
    ax.bar(x + width / 2, df["n_sig_cluster"], width, label="cluster", color="crimson")
    for xi, (_, row) in zip(x, df.iterrows()):
        ax.text(xi - width / 2, row["n_sig_supercluster"] + 0.15, str(int(row["n_sig_supercluster"])), ha="center", fontsize=8)
        ax.text(xi + width / 2, row["n_sig_cluster"] + 0.15, str(int(row["n_sig_cluster"])), ha="center", fontsize=8)
    ax.set_xticks(x, traits)
    ax.set_ylabel(f"Programs with FDR q < {FDR_ALPHA}")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, out_path, title)


def plot_cluster_forest(
    df: pd.DataFrame,
    out_path: Path,
    title: str,
    max_rows: int = 12,
) -> Path:
    """Forest plot of cluster-level coefficients with 95% confidence intervals.

    FDR-significant programs are preferred for display; if a table has no
    FDR hits, its strongest nominal programs are shown in grey instead. This
    function plots committed estimates only and recomputes no inference.

    Args:
        df: Table from :func:`load_enrichment_table`.
        out_path: Destination figure path.
        title: Axes title, naming trait and level.
        max_rows: Maximum number of programs to display.

    Returns:
        The written path.
    """
    hits = df[df["q_value"] < FDR_ALPHA]
    display = (hits if len(hits) else df).head(max_rows).iloc[::-1]
    y = np.arange(len(display))
    coefficients = display["coefficient"].to_numpy(dtype=float)
    ci95 = 1.96 * display["standard_error"].to_numpy(dtype=float)
    colors = np.where(display["q_value"] < FDR_ALPHA, "crimson", "0.6")

    fig_height = max(2.6, 0.34 * len(display) + 1.2)
    fig, ax = plt.subplots(figsize=(6.4, fig_height))
    ax.errorbar(
        coefficients,
        y,
        xerr=ci95,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=0.9,
        capsize=2,
        markersize=3.5,
    )
    ax.scatter(coefficients, y, c=colors, s=20, zorder=3)
    ax.axvline(0.0, color="0.3", linewidth=0.6)
    ax.set_yticks(y, display["program"].astype(str))
    ax.set_xlabel("Program coefficient (beta, 95% CI)")
    ax.set_title(title)
    ax.margins(x=0.18)
    for yi, (_, row) in zip(y, display.iterrows()):
        ax.text(
            1.01,
            yi,
            f"q={row['q_value']:.1e}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=7,
        )
    return _save(fig, out_path, title)


def discover_traits(reports_dir: Path) -> list[str]:
    """Trait tags of committed enrichment tables, in stable sorted order."""
    tags = set()
    for path in Path(reports_dir).glob("magma_enrichment_*.csv"):
        m = re.fullmatch(r"magma_enrichment_(.+)_(class|supercluster|cluster)\.csv", path.name)
        if m:
            tags.add(m.group(1))
    return sorted(tags)


def build_atlas_figures(
    reports_dir: Path,
    out_dir: Path,
    formats: tuple[str, ...] = FIGURE_FORMATS,
) -> list[Path]:
    """Regenerate every enrichment atlas figure from committed tables.

    Missing (trait x level) tables are skipped with a warning rather than
    treated as an error, so the known PD cluster-level gap (issue #10) does
    not block the committed cells.

    Args:
        reports_dir: Directory of committed ``reports/`` tables.
        out_dir: Figure destination (conventionally ``figures/``).
        formats: Output formats to write. Supported values are ``png`` and
            ``svg``; the default writes both.

    Returns:
        Sorted list of written figure paths.

    Raises:
        FileNotFoundError: If ``reports_dir`` or ``resolution_gain.csv``
            does not exist.
        ValueError: If an unsupported output format is requested.
    """
    reports_dir = Path(reports_dir)
    out_dir = Path(out_dir)
    formats = tuple(dict.fromkeys(formats))
    unsupported = set(formats) - set(FIGURE_FORMATS)
    if unsupported:
        raise ValueError(f"unsupported figure formats: {sorted(unsupported)}")
    if not reports_dir.is_dir():
        raise FileNotFoundError(f"reports directory not found: {reports_dir}")
    gain_path = reports_dir / "resolution_gain.csv"
    if not gain_path.is_file():
        raise FileNotFoundError(f"resolution-gain table not found: {gain_path}")

    written: list[Path] = []
    traits = discover_traits(reports_dir)
    for trait_index, trait in enumerate(traits):
        for level in LEVELS:
            path = enrichment_table_path(reports_dir, trait, level)
            if not path.is_file():
                logger.warning("skipping %s x %s: %s not committed", trait, level, path.name)
                continue
            df = load_enrichment_table(path)
            title = f"{trait} - {level} level"
            written.extend(
                _write_formats(
                    plot_volcano,
                    df,
                    base_path=out_dir / f"volcano_{trait}_{level}",
                    formats=formats,
                    title=title,
                )
            )
            written.extend(
                _write_formats(
                    plot_top_programs,
                    df,
                    base_path=out_dir / f"top_programs_{trait}_{level}",
                    formats=formats,
                    title=title,
                )
            )
            if level == "cluster":
                figure_number = 3 + trait_index
                written.extend(
                    _write_formats(
                        plot_cluster_forest,
                        df,
                        base_path=out_dir / f"figure{figure_number}_{trait}_cluster_forest",
                        formats=formats,
                        svg_plotter=plot_cluster_forest_svg,
                        title=f"{trait} cluster-level enrichment",
                    )
                )
    gain = load_resolution_gain(gain_path)
    written.extend(
        _write_formats(
            plot_atlas_overview,
            gain,
            base_path=out_dir / "figure1_atlas_overview",
            formats=formats,
            svg_plotter=plot_atlas_overview_svg,
        )
    )
    written.extend(
        _write_formats(
            plot_resolution_gain,
            gain,
            base_path=out_dir / "figure2_resolution_gain",
            formats=formats,
            svg_plotter=plot_resolution_gain_svg,
        )
    )
    # Backwards-compatible name used by the 2026-11-20 verification log.
    written.extend(
        _write_formats(
            plot_resolution_gain,
            gain,
            base_path=out_dir / "resolution_gain",
            formats=formats,
        )
    )
    return sorted(written)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: regenerate all atlas figures from ``reports/``."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reports-dir", type=Path, default=Path("reports"))
    parser.add_argument("--out-dir", type=Path, default=Path("figures"))
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=FIGURE_FORMATS,
        default=list(FIGURE_FORMATS),
        help="figure formats to write (default: png svg)",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    written = build_atlas_figures(args.reports_dir, args.out_dir, tuple(args.formats))
    for path in written:
        logger.info("wrote %s", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
