import numpy as np
import pandas as pd
import pytest

from celltype_heritability.magma import (
    bh_fdr,
    gene_test_statistics,
    map_snps_to_genes,
    program_enrichment,
)
from celltype_heritability.simulate import simulate_gwas


def test_magma_recovers_planted_enrichment():
    ds = simulate_gwas(seed=2)
    snp_map = map_snps_to_genes(ds.snps, ds.genes, window=2_000)
    assert set(snp_map) == set(ds.genes["gene"])
    gene_stats = gene_test_statistics(ds.z, snp_map)
    sizes = pd.Series(
        {g: len(snp_map[g]) for g in gene_stats.index}, index=gene_stats.index
    )
    res = program_enrichment(gene_stats, ds.gene_programs, gene_sizes=sizes)
    assert res["coefficient"].idxmax() == ds.enriched_celltype
    assert res.loc[ds.enriched_celltype, "p_value"] < 0.01


def test_program_enrichment_works_without_gene_sizes():
    """The documented ``gene_sizes=None`` default must not raise.

    A constant ``gene_sizes`` makes ``log(n_snps)`` an all-zero column, which
    is collinear with the intercept and used to raise ``LinAlgError``.
    """
    ds = simulate_gwas(seed=2)
    snp_map = map_snps_to_genes(ds.snps, ds.genes, window=2_000)
    gene_stats = gene_test_statistics(ds.z, snp_map)
    res = program_enrichment(gene_stats, ds.gene_programs)
    assert set(res.index) == set(ds.gene_programs)
    assert res["coefficient"].notna().all()
    # The planted program must still be the strongest one without the covariate.
    assert res["coefficient"].idxmax() == ds.enriched_celltype


def test_program_enrichment_reports_untestable_programs_as_nan():
    """Degenerate programs are unidentifiable, so they yield NaN, not a crash.

    A program with no tested genes (all-zero indicator) and a program covering
    every tested gene (all-one indicator) are both collinear with the
    intercept, so the model cannot separate them.
    """
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(40)]
    gene_stats = pd.Series(rng.chisquare(1.0, 40), index=genes)
    res = program_enrichment(
        gene_stats,
        {
            "empty": ["not_a_gene_1", "not_a_gene_2"],
            "everything": genes,
            "real": genes[:8],
        },
    )
    assert res.loc["empty", "coefficient"] != res.loc["empty", "coefficient"]
    assert res.loc["empty", "coefficient"] != res.loc["empty", "coefficient"]
    assert res.loc["everything", "coefficient"] != res.loc["everything", "coefficient"]
    assert np.isfinite(res.loc["real", "coefficient"])


def test_bh_fdr_ignores_nan_p_values():
    """One untestable program must not erase every other q-value.

    ``nan`` sorts last in numpy, so a naive BH step propagates NaN through the
    reverse cumulative minimum and blanks the whole family.
    """
    p = pd.Series({"a": 0.001, "b": 0.01, "c": 0.02, "untestable": np.nan})
    q = bh_fdr(p)
    assert q["untestable"] != q["untestable"]  # NaN in, NaN out
    for name in ("a", "b", "c"):
        assert np.isfinite(q[name]), f"q-value for {name} was poisoned by a NaN"
    # The three testable p-values form their own family of size 3:
    # 0.001*3/1, 0.01*3/2, 0.02*3/3, already monotone increasing.
    assert q["a"] == pytest.approx(0.003)
    assert q["b"] == pytest.approx(0.015)
    assert q["c"] == pytest.approx(0.02)


def test_bh_fdr_matches_hand_computed_values():
    """BH on a complete, NaN-free family is unchanged by the NaN handling."""
    p = pd.Series([0.01, 0.02, 0.03, 0.04, 0.05])
    q = bh_fdr(p)
    ranked = np.sort(p.to_numpy())
    n = len(ranked)
    expected = ranked * n / np.arange(1, n + 1)
    expected = np.minimum.accumulate(expected[::-1])[::-1]
    assert np.allclose(q.to_numpy(), np.minimum(expected, 1.0))
