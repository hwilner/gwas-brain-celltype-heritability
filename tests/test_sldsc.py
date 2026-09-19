import numpy as np

from celltype_heritability.ldscore import compute_ld_scores
from celltype_heritability.simulate import simulate_gwas
from celltype_heritability.sldsc import run_sldsc_regression


def test_sldsc_recovers_planted_enrichment():
    ds = simulate_gwas(seed=1)
    l2 = compute_ld_scores(ds.genotypes, ds.snps["pos"].to_numpy(), ds.annotations)
    res = run_sldsc_regression(ds.z**2, l2, n_samples=50_000.0)
    taus = res.coefficients.drop("intercept")
    # planted category has the largest positive coefficient
    assert taus.idxmax() == ds.enriched_celltype
    assert taus[ds.enriched_celltype] > 0
    assert res.z_scores[ds.enriched_celltype] > 2.0
    # non-planted categories are not enriched (one-sided check)
    others = res.z_scores.drop(["intercept", ds.enriched_celltype])
    assert (others < 2.0).all()
