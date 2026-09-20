#!/usr/bin/env python3
"""Run the real MAGMA-style cell-type enrichment atlas (issue #1).

Inputs (staged per docs/DATA_SOURCES.md; raw files are NOT committed):

* ``--atlas-dir`` : ABC Atlas tables (see ``celltype_heritability.atlas``)
* ``--gwas``      : name=path pairs of raw GWAS TSV(.gz) from GWAS Catalog
* ``--genes``     : Ensembl GRCh38 BioMart gene table
* ``--out-dir``   : reports/ output directory

Per trait: stream-harmonize summary statistics, map SNPs to genes (10 kb
flank), compute gene-level mean chi-square statistics, regress on each
marker program, and BH-FDR correct within each trait x level family. Also
runs the top_n=200 sensitivity for supercluster concordance and quantifies
the cluster-vs-supercluster resolution gain. Writes small CSV/JSON tables.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from celltype_heritability import atlas, magma, sumstats


def _parse_gwas(pairs: list[str]) -> dict[str, str]:
    return dict(p.split("=", 1) for p in pairs)


def enrichment_table(
    df: pd.DataFrame,
    programs: dict[str, list[str]],
    genes: pd.DataFrame,
    window: int = 10_000,
) -> pd.DataFrame:
    snps = df[["chrom", "pos"]].copy()
    snp_map = magma.map_snps_to_genes(snps, genes, window=window)
    gene_stats = magma.gene_test_statistics(df["z"].to_numpy(), snp_map)
    sizes = pd.Series(
        {g: len(snp_map[g]) for g in gene_stats.index}, index=gene_stats.index
    )
    progs = {k: v for k, v in programs.items() if len(v) >= 10}
    res = magma.program_enrichment(gene_stats, progs, gene_sizes=sizes)
    res.index.name = "program"
    res = res.reset_index()
    res["q_value"] = magma.bh_fdr(res["p_value"])
    res["n_program_genes"] = [len(programs[p]) for p in res["program"]]
    return res.sort_values("p_value").reset_index(drop=True)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--atlas-dir", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--gwas", nargs="+", required=True, help="NAME=path.tsv.gz")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=100)
    args = p.parse_args(argv)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genes = sumstats.load_gene_annotation(args.genes)
    gwas = _parse_gwas(args.gwas)

    # ---- gene programs at three hierarchy levels (committed, small)
    levels = ["class", "supercluster", "cluster"]
    manifest_md5 = json.loads(
        (Path(args.atlas_dir) / "manifest.json").read_text()
    )
    programs = {}
    for level in levels:
        programs[level] = atlas.build_marker_programs(
            args.atlas_dir, level=level, top_n=args.top_n
        )
        atlas.export_programs(
            programs[level], out / "gene_programs", manifest_md5, level, args.top_n
        )

    meta: dict = {
        "pipeline": "in-repo MAGMA-style (scripts/run_magma_enrichment.py)",
        "window_bp": 10_000,
        "top_n": args.top_n,
        "n_genes": int(len(genes)),
        "atlas_md5": manifest_md5,
        "traits": {},
    }

    results = {}
    harmonized = {}
    for trait, path in gwas.items():
        src = sumstats.GWAS_SOURCES.get(trait, {})
        df, qc = sumstats.stream_harmonize(
            path,
            n_cases=float(src["n_cases"]) if src.get("n_cases") else None,
            n_controls=float(src["n_controls"]) if src.get("n_controls") else None,
            minimal=True,
        )
        harmonized[trait] = df
        meta["traits"].setdefault(trait, {})["qc"] = qc
        meta["traits"][trait]["gwas_md5"] = sumstats.md5_of(path)
    for trait in gwas:
        df = harmonized[trait]
        for level in levels:
            res = enrichment_table(df, programs[level], genes)
            res["trait"] = trait
            res["level"] = level
            results[(trait, level)] = res
            res.to_csv(out / f"magma_enrichment_{trait}_{level}.csv", index=False)
            meta["traits"][trait][level] = {
                "n_programs_tested": int(len(res)),
                "n_fdr05": int((res["q_value"] < 0.05).sum()),
                "top": res.iloc[0]["program"] if len(res) else None,
            }

    # ---- concordance: supercluster top_n=200 sensitivity
    sens_rows = []
    for trait in gwas:
        progs200 = atlas.build_marker_programs(
            args.atlas_dir, level="supercluster", top_n=200
        )
        res200 = enrichment_table(harmonized[trait], progs200, genes)
        res100 = results[(trait, "supercluster")]
        merged = res100.merge(
            res200[["program", "z_score"]], on="program", suffixes=("_100", "_200")
        )
        rho = stats.spearmanr(merged["z_score_100"], merged["z_score_200"]).statistic
        sens_rows.append({"trait": trait, "spearman_z_top100_vs_top200": float(rho)})
    pd.DataFrame(sens_rows).to_csv(out / "concordance_topn_sensitivity.csv", index=False)

    # ---- resolution gain: cluster vs supercluster
    gain = []
    for trait in gwas:
        c = results[(trait, "cluster")]
        s = results[(trait, "supercluster")]
        gain.append(
            {
                "trait": trait,
                "n_sig_supercluster": int((s["q_value"] < 0.05).sum()),
                "n_sig_cluster": int((c["q_value"] < 0.05).sum()),
                "min_p_supercluster": float(s["p_value"].min()),
                "min_p_cluster": float(c["p_value"].min()),
            }
        )
    pd.DataFrame(gain).to_csv(out / "resolution_gain.csv", index=False)

    (out / "RUN_METADATA.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps({t: meta["traits"][t] for t in meta["traits"]}, indent=2))


if __name__ == "__main__":
    main()
