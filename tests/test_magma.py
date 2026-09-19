import numpy as np
import pandas as pd

from celltype_heritability.magma import (
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
