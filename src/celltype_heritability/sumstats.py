"""Download, harmonize, and QC public GWAS summary statistics.

Staged datasets (all public, GRCh38):

* AD  - Bellenguez et al. 2022 (Nat Genet), EBI GWAS Catalog GCST90027158
        (GRCh38-native, 21,101,114 variants) and FinnGen R11 G6_ALZHEIMER.
* PD  - Nalls et al. 2019 (Lancet Neurol) excluding 23andMe, EBI GWAS
        Catalog GCST009325, harmonised build (GRCh38), and FinnGen R11
        G6_PARKINSON.

The GWAS Catalog files are large (0.5-0.8 GB gzipped); ``download``/
``stream_harmonize`` read them in chunks so they can be processed line-by-line
without loading the raw file into memory. Only the small harmonized/QC'd
outputs are kept for analysis.

A small Ensembl (GRCh38) gene annotation, fetched via BioMart, provides the
SNP -> gene map used by the MAGMA-style analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

#: Registry of staged GWAS summary statistics.
GWAS_SOURCES: dict[str, dict[str, str]] = {
    "AD_Bellenguez2022": {
        "accession": "GCST90027158",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz",
        "build": "GRCh38",
        "citation": "Bellenguez et al. 2022, Nat Genet (PMID 35379992)",
    },
    "PD_Nalls2019": {
        "accession": "GCST009325",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz",
        "build": "GRCh38 (harmonised)",
        "citation": "Nalls et al. 2019, Lancet Neurol (PMID 31501809)",
    },
    "AD_FinnGen_R11": {
        "accession": "finngen_R11_G6_ALZHEIMER",
        "url": "https://storage.googleapis.com/finngen-public-data-r11/"
        "summary_stats/finngen_R11_G6_ALZHEIMER.gz",
        "build": "GRCh38",
        "citation": "FinnGen release 11 (2024), endpoint G6_ALZHEIMER; "
        "11,755 cases / 441,978 controls",
        "n_cases": "11755",
        "n_controls": "441978",
    },
    "PD_FinnGen_R11": {
        "accession": "finngen_R11_G6_PARKINSON",
        "url": "https://storage.googleapis.com/finngen-public-data-r11/"
        "summary_stats/finngen_R11_G6_PARKINSON.gz",
        "build": "GRCh38",
        "citation": "FinnGen release 11 (2024), endpoint G6_PARKINSON; "
        "5,150 cases / 448,583 controls",
        "n_cases": "5150",
        "n_controls": "448583",
    },
}

#: Ensembl BioMart query (GRCh38) for the gene annotation used genome-wide.
ENSEMBL_BIOMART_URL = "https://useast.ensembl.org/biomart/martservice"
ENSEMBL_GENE_QUERY = (
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
    '<Query virtualSchemaName="default" formatter="TSV" header="1" '
    'uniqueRows="0" datasetConfigVersion="0.6">'
    '<Dataset name="hsapiens_gene_ensembl" interface="default">'
    '<Attribute name="ensembl_gene_id"/>'
    '<Attribute name="external_gene_name"/>'
    '<Attribute name="chromosome_name"/>'
    '<Attribute name="start_position"/>'
    '<Attribute name="end_position"/>'
    '<Attribute name="gene_biotype"/>'
    "</Dataset></Query>"
)

CHUNK = 1 << 20
AUTOSOMES = {str(c) for c in range(1, 23)}

# Flexible column look-ups across raw/harmonised GWAS Catalog formats.
_COLS = {
    "rsid": ["rsid", "variant_id", "hm_rsid", "SNP", "rsids"],
    "chrom": ["chromosome", "hm_chrom", "chrom", "#chrom"],
    "pos": ["base_pair_location", "hm_pos", "pos"],
    "effect_allele": ["effect_allele", "hm_effect_allele", "alt"],
    "other_allele": ["other_allele", "hm_other_allele", "ref"],
    "beta": ["beta", "hm_beta"],
    "se": ["standard_error", "hm_standard_error", "sebeta"],
    "p": ["p_value", "hm_p_value", "pval"],
    "n_cases": ["n_cases", "N_cases"],
    "n_controls": ["n_controls", "N_controls"],
}


def download_gene_annotation(dest: str | Path) -> pd.DataFrame:
    """Fetch the Ensembl GRCh38 gene table via BioMart (POST) and save TSV."""
    import urllib.parse
    import urllib.request

    req = urllib.request.Request(
        ENSEMBL_BIOMART_URL,
        data=urllib.parse.urlencode({"query": ENSEMBL_GENE_QUERY}).encode(),
    )
    with urllib.request.urlopen(req) as resp:
        data = resp.read()
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return load_gene_annotation(dest)


def load_gene_annotation(path: str | Path) -> pd.DataFrame:
    """Load the BioMart gene table -> DataFrame(gene, chrom, start, end).

    Protein-coding autosomal genes, longest span kept per gene symbol.
    """
    raw = pd.read_csv(path, sep="\t").rename(
        columns={
            "Gene stable ID": "ensembl",
            "Gene name": "gene",
            "Chromosome/scaffold name": "chrom",
            "Gene start (bp)": "start",
            "Gene end (bp)": "end",
            "Gene type": "biotype",
        }
    )
    df = raw[raw["biotype"] == "protein_coding"].copy()
    df = df[df["chrom"].astype(str).isin(AUTOSOMES)]
    df = df.dropna(subset=["gene"])
    df["chrom"] = df["chrom"].astype(int)
    df["span"] = df["end"] - df["start"]
    df = df.sort_values("span", ascending=False).drop_duplicates("gene")
    return (
        df[["gene", "chrom", "start", "end"]]
        .sort_values(["chrom", "start"])
        .reset_index(drop=True)
    )


def _resolve_columns(header: list[str]) -> dict[str, str]:
    lut = {}
    lower = {c.lower(): c for c in header}
    for std, candidates in _COLS.items():
        for cand in candidates:
            if cand.lower() in lower:
                lut[std] = lower[cand.lower()]
                break
    return lut


def stream_harmonize(
    path: str | Path,
    chunksize: int = 500_000,
    max_rows: int | None = None,
    n_cases: float | None = None,
    n_controls: float | None = None,
    minimal: bool = False,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Stream a GWAS Catalog TSV(.gz) and harmonize to a common schema.

    Output columns: ``rsid, chrom, pos, effect_allele, other_allele, beta,
    se, p, z, n``. With ``minimal=True`` only ``chrom, pos, z`` are kept
    (downcast dtypes) for genome-wide files in memory-tight environments.
    QC (applied row-wise, deterministic): autosomes 1-22, valid p in (0, 1],
    non-missing beta/se with se > 0, standard ACGT SNP alleles. Returns
    (harmonized frame, QC counts).
    """
    counts = {"rows_in": 0, "autosomes": 0, "valid_stats": 0}
    out_chunks = []
    reader = pd.read_csv(
        path, sep="\t", chunksize=chunksize, dtype=str, low_memory=False
    )
    lut = None
    for chunk in reader:
        if lut is None:
            lut = _resolve_columns(list(chunk.columns))
            missing = {"chrom", "pos", "beta", "se", "p"} - set(lut)
            if missing:
                raise ValueError(f"cannot resolve required columns {sorted(missing)}")
        counts["rows_in"] += len(chunk)
        c = chunk.rename(columns={v: k for k, v in lut.items()})
        c = c[c["chrom"].isin(AUTOSOMES)]
        counts["autosomes"] += len(c)
        for col in ("pos", "beta", "se", "p"):
            c[col] = pd.to_numeric(c[col], errors="coerce")
        keep = (
            c["pos"].notna()
            & c["beta"].notna()
            & c["se"].notna()
            & (c["se"] > 0)
            & c["p"].notna()
            & (c["p"] > 0)
            & (c["p"] <= 1)
        )
        if "effect_allele" in c and "other_allele" in c:
            keep &= c["effect_allele"].isin(list("ACGT")) & c["other_allele"].isin(
                list("ACGT")
            )
        c = c[keep]
        counts["valid_stats"] += len(c)
        if not len(c):
            continue
        if {"n_cases", "n_controls"} <= set(lut):
            n = pd.to_numeric(c.get("n_cases"), errors="coerce") + pd.to_numeric(
                c.get("n_controls"), errors="coerce"
            )
        elif n_cases is not None and n_controls is not None:
            n = pd.Series(float(n_cases) + float(n_controls), index=c.index)
        else:
            n = pd.Series(np.nan, index=c.index)
        if minimal:
            out_chunks.append(
                pd.DataFrame(
                    {
                        "chrom": c["chrom"].astype("int8"),
                        "pos": c["pos"].astype("int32"),
                        "z": (c["beta"] / c["se"]).astype("float32"),
                    }
                )
            )
        else:
            out_chunks.append(
                pd.DataFrame(
                    {
                        "rsid": c.get("rsid", pd.Series(".", index=c.index)),
                        "chrom": c["chrom"],
                        "pos": c["pos"].astype(int),
                        "effect_allele": c.get(
                            "effect_allele", pd.Series(".", index=c.index)
                        ),
                        "other_allele": c.get(
                            "other_allele", pd.Series(".", index=c.index)
                        ),
                        "beta": c["beta"],
                        "se": c["se"],
                        "p": c["p"],
                        "z": c["beta"] / c["se"],
                        "n": n,
                    }
                )
            )
        if max_rows is not None and counts["rows_in"] >= max_rows:
            break
    df = (
        pd.concat(out_chunks, ignore_index=True)
        if out_chunks
        else pd.DataFrame(
            columns=[
                "rsid", "chrom", "pos", "effect_allele", "other_allele",
                "beta", "se", "p", "z", "n",
            ]
        )
    )
    before = len(df)
    df = df.drop_duplicates(subset=["chrom", "pos"], keep="first").reset_index(drop=True)
    counts["dedup_dropped"] = before - len(df)
    counts["rows_out"] = len(df)
    return df, counts


def md5_of(path: str | Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("download", help="download a registered GWAS (streamed)")
    d.add_argument("name", choices=sorted(GWAS_SOURCES))
    d.add_argument("--dest", required=True)
    h = sub.add_parser("harmonize", help="harmonize a raw GWAS TSV(.gz)")
    h.add_argument("--src", required=True)
    h.add_argument("--dest", required=True)
    h.add_argument("--max-rows", type=int, default=None)
    args = p.parse_args(argv)
    if args.cmd == "download":
        from .atlas import download_file

        md5 = download_file(GWAS_SOURCES[args.name]["url"], args.dest)
        print(json.dumps({"file": args.dest, "md5": md5}))
    else:
        df, counts = stream_harmonize(args.src, max_rows=args.max_rows)
        df.to_csv(args.dest, sep="\t", index=False, compression="gzip")
        print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
