"""Integrity checks for the committed Paper 1 results package (issue #8).

These tests are deliberately cheap and data-free: they verify that the small
committed result tables agree with the pinned run metadata and that the
committed reviewable figure files stay wired to their captions. They do not
recompute any scientific result and they do not replace the external-data
clean-room re-run still required by issue #8.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
LEVELS = ("class", "supercluster", "cluster")
FDR_ALPHA = 0.05


def _metadata() -> dict:
    return json.loads((ROOT / "reports" / "RUN_METADATA.json").read_text())


def _table_path(trait: str, level: str) -> Path:
    return REPORTS / f"magma_enrichment_{trait}_{level}.csv"


def _assert_table_matches_metadata(trait: str, level: str) -> None:
    """Check one committed table against its pinned RUN_METADATA entry."""
    path = _table_path(trait, level)
    assert path.is_file(), f"expected committed table is missing: {path.name}"
    table = pd.read_csv(path)
    expected = _metadata()["traits"][trait][level]
    assert len(table) == expected["n_programs_tested"]
    assert int((table["q_value"] < FDR_ALPHA).sum()) == expected["n_fdr05"]
    assert table.iloc[0]["program"] == expected["top"]
    assert table["p_value"].is_monotonic_increasing


def test_committed_enrichment_tables_match_run_metadata():
    """Every committed trait x level table agrees with RUN_METADATA."""
    metadata = _metadata()
    committed = 0
    for trait in metadata["traits"]:
        for level in LEVELS:
            if _table_path(trait, level).is_file():
                _assert_table_matches_metadata(trait, level)
                committed += 1
    # Five tables are committed at the time of this test; the PD cluster table
    # is the single known gap tracked by issue #10. Once it lands this count
    # should become six without any other test change.
    assert committed >= 5


def test_pd_cluster_gap_is_explicit_until_the_table_lands():
    """The one uncommitted table must stay explicitly marked, not silent."""
    metadata = _metadata()
    expected = metadata["traits"]["PD_FinnGen_R11"]["cluster"]
    results_text = (ROOT / "docs" / "RESULTS.md").read_text()
    path = _table_path("PD_FinnGen_R11", "cluster")
    if path.is_file():
        _assert_table_matches_metadata("PD_FinnGen_R11", "cluster")
        assert "the PD cluster-level full table\n   (`magma_enrichment_PD_FinnGen_R11_cluster.csv`) is not committed" not in results_text
    else:
        assert expected == {
            "n_programs_tested": 436,
            "n_fdr05": 11,
            "top": "CA4_197",
        }
        assert "Known gap" in results_text
        assert "issue #10" in results_text


def test_resolution_gain_matches_committed_tables_and_metadata():
    """The headline resolution table agrees with its committed sources."""
    metadata = _metadata()
    gain = pd.read_csv(REPORTS / "resolution_gain.csv").set_index("trait")
    assert set(gain.index) == set(metadata["traits"])
    for trait, row in gain.iterrows():
        supercluster = pd.read_csv(_table_path(trait, "supercluster"))
        assert row["n_sig_supercluster"] == int(
            (supercluster["q_value"] < FDR_ALPHA).sum()
        )
        assert row["min_p_supercluster"] == pytest.approx(
            supercluster["p_value"].min()
        )
        cluster_path = _table_path(trait, "cluster")
        if cluster_path.is_file():
            cluster = pd.read_csv(cluster_path)
            assert row["n_sig_cluster"] == int((cluster["q_value"] < FDR_ALPHA).sum())
            assert row["min_p_cluster"] == pytest.approx(cluster["p_value"].min())
        else:
            # The PD cluster cell is still sourced from pinned run metadata
            # until issue #10 commits the full table.
            expected = metadata["traits"][trait]["cluster"]
            assert row["n_sig_cluster"] == expected["n_fdr05"]


def test_committed_publication_svgs_exist_and_are_captioned():
    """The text-format publication figures stay committed and discoverable."""
    results_text = (ROOT / "docs" / "RESULTS.md").read_text()
    expected = [
        "figure1_atlas_overview.svg",
        "figure2_resolution_gain.svg",
        "figure3_AD_FinnGen_R11_cluster_forest.svg",
    ]
    pd_cluster = _table_path("PD_FinnGen_R11", "cluster")
    if pd_cluster.is_file():
        expected.append("figure4_PD_FinnGen_R11_cluster_forest.svg")
    for name in expected:
        path = ROOT / "figures" / name
        assert path.is_file(), f"committed publication figure missing: {name}"
        text = path.read_text()
        assert "<svg" in text[:500]
        assert f"figures/{name}" in results_text
