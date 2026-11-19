"""Tests for the enrichment atlas figure scripts (issue #7).

Fully synthetic and data-free per CONTRIBUTING.md: a tiny enrichment table
is fabricated in ``tmp_path`` and the tests assert that figures are written,
non-trivial in size, deterministic across runs, and that missing tables are
skipped rather than crashing the build.
"""

import logging

import pandas as pd
import pytest

from celltype_heritability.figures import (
    build_atlas_figures,
    discover_traits,
    enrichment_table_path,
    load_enrichment_table,
)


def _synthetic_table(programs=("Astro", "Mgl", "OPC", "Vascular")) -> pd.DataFrame:
    """Four-program synthetic enrichment table with one planted FDR hit."""
    return pd.DataFrame(
        {
            "program": list(programs),
            "coefficient": [1.5, 0.4, -0.2, 0.9],
            "standard_error": [0.2, 0.25, 0.22, 0.3],
            "z_score": [7.5, 1.6, -0.9, 3.0],
            "p_value": [1e-12, 0.11, 0.36, 0.0027],
            "n_genes": [100, 100, 100, 100],
            "q_value": [4e-12, 0.18, 0.36, 0.0054],
            "n_program_genes": [100, 100, 100, 100],
            "trait": ["SYNTH"] * 4,
            "level": ["supercluster"] * 4,
        }
    )


def _write_reports(reports_dir) -> None:
    """Stage a minimal synthetic reports/ directory."""
    reports_dir.mkdir()
    _synthetic_table().to_csv(
        enrichment_table_path(reports_dir, "SYNTH", "supercluster"), index=False
    )
    pd.DataFrame(
        {
            "trait": ["SYNTH"],
            "n_sig_supercluster": [1],
            "n_sig_cluster": [3],
            "min_p_supercluster": [1e-12],
            "min_p_cluster": [1e-20],
        }
    ).to_csv(reports_dir / "resolution_gain.csv", index=False)


def test_load_enrichment_table_validates_schema(tmp_path):
    bad = tmp_path / "bad.csv"
    _synthetic_table().drop(columns=["q_value"]).to_csv(bad, index=False)
    with pytest.raises(ValueError, match="q_value"):
        load_enrichment_table(bad)


def test_build_atlas_figures_writes_expected_pngs(tmp_path):
    reports = tmp_path / "reports"
    _write_reports(reports)
    out = tmp_path / "figures"
    written = build_atlas_figures(reports, out)
    names = {p.name for p in written}
    assert names == {
        "volcano_SYNTH_supercluster.png",
        "top_programs_SYNTH_supercluster.png",
        "resolution_gain.png",
    }
    for path in written:
        # A real PNG with content, not an empty placeholder.
        assert path.stat().st_size > 5_000


def test_build_atlas_figures_is_deterministic(tmp_path):
    reports = tmp_path / "reports"
    _write_reports(reports)
    first = build_atlas_figures(reports, tmp_path / "run1")
    second = build_atlas_figures(reports, tmp_path / "run2")
    assert [p.name for p in first] == [p.name for p in second]
    for a, b in zip(first, second):
        assert a.read_bytes() == b.read_bytes(), f"{a.name} is not deterministic"


def test_missing_tables_are_skipped_with_warning(tmp_path, caplog):
    """The PD cluster-level gap (issue #10) must not abort the build."""
    reports = tmp_path / "reports"
    _write_reports(reports)
    with caplog.at_level(logging.WARNING):
        written = build_atlas_figures(reports, tmp_path / "figures")
    assert any("not committed" in r.message for r in caplog.records)
    assert {p.name for p in written} == {
        "volcano_SYNTH_supercluster.png",
        "top_programs_SYNTH_supercluster.png",
        "resolution_gain.png",
    }


def test_discover_traits_ignores_unrelated_csvs(tmp_path):
    reports = tmp_path / "reports"
    _write_reports(reports)
    (reports / "unrelated.csv").write_text("a,b\n1,2\n")
    assert discover_traits(reports) == ["SYNTH"]
