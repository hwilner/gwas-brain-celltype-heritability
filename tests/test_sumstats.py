import pandas as pd
import pytest

from celltype_heritability import sumstats


def _write_gwas(tmp_path, rows):
    p = tmp_path / "gwas.tsv.gz"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False, compression="gzip")
    return p


def _row(**kw):
    base = {
        "variant_id": "rs1",
        "p_value": 0.05,
        "chromosome": "1",
        "base_pair_location": 1000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.1,
        "standard_error": 0.05,
        "n_cases": 100,
        "n_controls": 1000,
    }
    base.update(kw)
    return base


def test_stream_harmonize_qc_and_schema(tmp_path):
    rows = [
        _row(),  # keep
        _row(variant_id="rs2", base_pair_location=2000, chromosome="X"),  # non-autosome
        _row(variant_id="rs3", base_pair_location=3000, p_value=0.0),  # bad p
        _row(variant_id="rs4", base_pair_location=4000, standard_error=0.0),  # bad se
        _row(variant_id="rs5", base_pair_location=1000),  # duplicate position
    ]
    df, counts = sumstats.stream_harmonize(_write_gwas(tmp_path, rows))
    assert counts["rows_in"] == 5
    assert counts["autosomes"] == 4
    assert counts["valid_stats"] == 2
    assert counts["rows_out"] == 1
    assert list(df.columns) == [
        "rsid", "chrom", "pos", "effect_allele", "other_allele",
        "beta", "se", "p", "z", "n",
    ]
    assert df.loc[0, "z"] == pytest.approx(2.0)
    assert df.loc[0, "n"] == 1100
    # determinism
    df2, counts2 = sumstats.stream_harmonize(_write_gwas(tmp_path, rows))
    assert df.equals(df2) and counts == counts2


def test_load_gene_annotation_filters_and_dedups(tmp_path):
    p = tmp_path / "genes.tsv"
    pd.DataFrame(
        {
            "Gene stable ID": ["E1", "E2", "E3", "E4"],
            "Gene name": ["G1", "G1", "G2", "G3"],
            "Chromosome/scaffold name": ["1", "1", "2", "MT"],
            "Gene start (bp)": [100, 100, 500, 1],
            "Gene end (bp)": [200, 300, 600, 50],
            "Gene type": ["protein_coding"] * 3 + ["Mt_rRNA"],
        }
    ).to_csv(p, sep="\t", index=False)
    genes = sumstats.load_gene_annotation(p)
    assert list(genes["gene"]) == ["G1", "G2"]  # longest G1 span kept; MT dropped
    assert genes.loc[genes["gene"] == "G1", "end"].iloc[0] == 300


def test_md5_of_stable(tmp_path):
    p = tmp_path / "f.bin"
    p.write_bytes(b"hello")
    assert sumstats.md5_of(p) == "5d41402abc4b2a76b9719d911017c592"
