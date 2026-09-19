import numpy as np
import pandas as pd

from celltype_heritability.sclinker import disease_relevance_scores, score_programs


def test_score_programs_deterministic_and_specific():
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(
        rng.normal(size=(10, 3)),
        index=[f"g{i}" for i in range(10)],
        columns=["astro", "neuron", "micro"],
    )
    expr.loc[["g0", "g1"], "micro"] += 5.0  # micro-specific markers
    programs = {"micro_prog": ["g0", "g1"], "other": ["g5", "g6"]}
    scores = score_programs(expr, programs)
    assert scores["micro_prog"].idxmax() == "micro"
    # determinism
    assert scores.equals(score_programs(expr, programs))


def test_disease_relevance_weights_by_gene_association():
    scores = pd.DataFrame(
        {"prog": [1.0, 2.0]}, index=["ct_a", "ct_b"]
    )
    gene_scores = pd.Series({"g1": 3.0}, dtype=float)
    out = disease_relevance_scores(scores, gene_scores, {"prog": ["g1"]})
    assert out.idxmax() == "ct_b"
    assert np.allclose(out, [3.0, 6.0])
