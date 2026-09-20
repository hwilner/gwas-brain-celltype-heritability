"""Allen Brain Cell Atlas (human whole-brain, Siletti et al. 2023) loaders.

Downloads the published ABC Atlas cell-type taxonomy and the MapMyCells
query-marker tables (ranked marker genes per taxonomy node, computed by the
Allen Institute ``cell_type_mapper`` v1.3.0, timestamp 2024-02-21) and builds
deterministic per-cell-type marker gene programs at three hierarchy levels:

* ``class``        : Neuronal vs Non-neuronal (from neurotransmitter calls)
* ``supercluster`` : the 31 WHB superclusters (e.g. "Microglia", "Astrocyte")
* ``cluster``      : the 461 WHB clusters (finer cell types)

Marker-selection rule (pre-registered, see docs/ANALYSIS_PLAN.md): for each
taxonomy node take the top ``top_n`` ranked MapMyCells query markers, map
Ensembl IDs to gene symbols with the WHB-10Xv3 gene table, keep protein-coding
genes only, and deduplicate preserving rank order. Class-level programs pool
cluster markers ranked by (frequency across member clusters, mean rank).

The full MapMyCells precomputed reference statistics
(``precomputed_stats.siletti.training.h5``, ~8.7 GB) are optional and only
fetched with ``full=True`` / the ``--full`` CLI flag; downloads are streamed
in chunks so large files never need to fit in memory at once.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ABC_S3 = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

#: Public ABC Atlas files used by the default (small-table) path.
ATLAS_FILES: dict[str, str] = {
    "query_markers.n10.20240221800.json": (
        f"{ABC_S3}/mapmycells/WHB-10Xv3/20240831/query_markers.n10.20240221800.json"
    ),
    "cluster.csv": f"{ABC_S3}/metadata/WHB-taxonomy/20240330/cluster.csv",
    "cluster_annotation_term.csv": (
        f"{ABC_S3}/metadata/WHB-taxonomy/20240330/cluster_annotation_term.csv"
    ),
    "cluster_annotation_term_set.csv": (
        f"{ABC_S3}/metadata/WHB-taxonomy/20240330/cluster_annotation_term_set.csv"
    ),
    "cluster_to_cluster_annotation_membership.csv": (
        f"{ABC_S3}/metadata/WHB-taxonomy/20240330/"
        "cluster_to_cluster_annotation_membership.csv"
    ),
    "gene.csv": f"{ABC_S3}/metadata/WHB-10Xv3/20241115/gene.csv",
}

#: Optional full MapMyCells reference statistics (~8.7 GB, streamed).
FULL_FILES: dict[str, str] = {
    "precomputed_stats.siletti.training.h5": (
        f"{ABC_S3}/mapmycells/WHB-10Xv3/20240831/"
        "precomputed_stats.siletti.training.h5"
    ),
}

LEVEL_TERM_SETS = {
    "supercluster": "CCN202210140_SUPC",
    "cluster": "CCN202210140_CLUS",
    "neurotransmitter": "CCN202210140_NEUR",
}

CHUNK = 1 << 20


def md5_of(path: str | Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, dest: str | Path) -> str:
    """Stream ``url`` to ``dest`` in 1 MB chunks; return the file's md5."""
    import urllib.request

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.md5()
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as out:
        while True:
            chunk = resp.read(CHUNK)
            if not chunk:
                break
            h.update(chunk)
            out.write(chunk)
    return h.hexdigest()


def download_atlas(data_dir: str | Path, full: bool = False) -> dict[str, str]:
    """Download the ABC Atlas tables (and optionally the full stats h5).

    Returns a mapping of file name -> md5 checksum (the run manifest).
    """
    files = dict(ATLAS_FILES)
    if full:
        files.update(FULL_FILES)
    manifest = {}
    for name, url in files.items():
        manifest[name] = download_file(url, Path(data_dir) / name)
    (Path(data_dir) / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def load_taxonomy(data_dir: str | Path) -> pd.DataFrame:
    """Load the WHB taxonomy terms with hierarchy level and parent label."""
    terms = pd.read_csv(Path(data_dir) / "cluster_annotation_term.csv")
    return terms


def load_query_markers(data_dir: str | Path) -> dict[str, list[str]]:
    """Load ranked MapMyCells query markers: term key -> Ensembl gene IDs."""
    with open(Path(data_dir) / "query_markers.n10.20240221800.json") as fh:
        raw = json.load(fh)
    return {k: v for k, v in raw.items() if "/" in k and isinstance(v, list)}


def load_gene_table(data_dir: str | Path) -> pd.DataFrame:
    """WHB-10Xv3 gene table (Ensembl ID -> symbol, biotype)."""
    return pd.read_csv(Path(data_dir) / "gene.csv")


def _ensembl_to_symbol(gene_table: pd.DataFrame) -> dict[str, str]:
    pc = gene_table[gene_table["biotype"] == "protein_coding"]
    return dict(zip(pc["gene_identifier"], pc["gene_symbol"]))


def build_marker_programs(
    data_dir: str | Path,
    level: str = "cluster",
    top_n: int = 100,
) -> dict[str, list[str]]:
    """Build marker gene programs (gene symbols) at one hierarchy level.

    Pre-registered rule: top ``top_n`` ranked MapMyCells query markers per
    taxonomy node, protein-coding only, rank order preserved. For
    ``level="class"`` (Neuronal/Non-neuronal) cluster markers are pooled and
    ranked by (frequency across member clusters, mean rank, symbol).
    """
    if level not in ("class", "supercluster", "cluster"):
        raise ValueError(f"unknown level {level!r}")
    markers = load_query_markers(data_dir)
    terms = load_taxonomy(data_dir)
    id2sym = _ensembl_to_symbol(load_gene_table(data_dir))

    def to_symbols(ids: list[str], n: int) -> list[str]:
        out: list[str] = []
        for g in ids:
            sym = id2sym.get(g)
            if sym and sym not in out:
                out.append(sym)
            if len(out) >= n:
                break
        return out

    if level in ("supercluster", "cluster"):
        term_set = LEVEL_TERM_SETS[level]
        names = dict(zip(terms["label"], terms["name"]))
        programs = {}
        for key, ids in markers.items():
            ts, label = key.split("/", 1)
            if ts != term_set:
                continue
            programs[str(names.get(label, label))] = to_symbols(ids, top_n)
        return dict(sorted(programs.items()))

    # class level: neuronal vs non-neuronal from neurotransmitter calls
    neuronal = _neuronal_cluster_labels(data_dir, terms)
    clus_level = LEVEL_TERM_SETS["cluster"]
    pools: dict[str, dict[str, list[int]]] = {"Neuronal": {}, "Non-neuronal": {}}
    for key, ids in markers.items():
        ts, label = key.split("/", 1)
        if ts != clus_level:
            continue
        klass = "Neuronal" if label in neuronal else "Non-neuronal"
        for rank, g in enumerate(ids):
            sym = id2sym.get(g)
            if not sym:
                continue
            pools[klass].setdefault(sym, []).append(rank)

    def pool_top(pool: dict[str, list[int]], n: int) -> list[str]:
        scored = sorted(
            pool.items(), key=lambda kv: (-len(kv[1]), sum(kv[1]) / len(kv[1]), kv[0])
        )
        return [sym for sym, _ in scored[:n]]

    return {klass: pool_top(pool, top_n) for klass, pool in pools.items()}


def _neuronal_cluster_labels(data_dir: str | Path, terms: pd.DataFrame) -> set[str]:
    """Cluster labels whose cells carry a real neurotransmitter assignment."""
    mem = pd.read_csv(
        Path(data_dir) / "cluster_to_cluster_annotation_membership.csv"
    )
    neur = mem[mem["cluster_annotation_term_set_label"] == LEVEL_TERM_SETS["neurotransmitter"]]
    neur_terms = terms[terms["cluster_annotation_term_set_label"] == LEVEL_TERM_SETS["neurotransmitter"]]
    none_label = neur_terms.loc[neur_terms["name"].isna(), "label"].iloc[0]
    neuronal_alias = set(
        neur.loc[neur["cluster_annotation_term_label"] != none_label, "cluster_alias"]
    )
    clus_mem = mem[mem["cluster_annotation_term_set_label"] == LEVEL_TERM_SETS["cluster"]]
    alias_to_clus = dict(zip(clus_mem["cluster_alias"], clus_mem["cluster_annotation_term_label"]))
    return {alias_to_clus[a] for a in neuronal_alias if a in alias_to_clus}


def export_programs(
    programs: dict[str, list[str]],
    out_dir: str | Path,
    source_md5: dict[str, str],
    level: str,
    top_n: int,
) -> Path:
    """Write programs JSON + a manifest CSV (one row per program)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "level": level,
        "top_n": top_n,
        "source_md5": source_md5,
        "programs": programs,
    }
    path = out_dir / f"programs_{level}_top{top_n}.json"
    path.write_text(json.dumps(payload, indent=2))
    manifest = pd.DataFrame(
        {
            "program": list(programs),
            "n_genes": [len(v) for v in programs.values()],
            "level": level,
        }
    )
    manifest.to_csv(out_dir / f"programs_{level}_top{top_n}_manifest.csv", index=False)
    return path


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("download", help="download ABC Atlas tables")
    d.add_argument("--data-dir", required=True)
    d.add_argument(
        "--full",
        action="store_true",
        help="also stream the ~8.7 GB MapMyCells precomputed stats h5",
    )
    b = sub.add_parser("build", help="build marker gene programs")
    b.add_argument("--data-dir", required=True)
    b.add_argument("--out-dir", required=True)
    b.add_argument("--level", default="cluster",
                   choices=["class", "supercluster", "cluster"])
    b.add_argument("--top-n", type=int, default=100)
    args = p.parse_args(argv)
    if args.cmd == "download":
        manifest = download_atlas(args.data_dir, full=args.full)
        for name, md5 in manifest.items():
            print(f"{md5}  {name}")
    else:
        manifest_path = Path(args.data_dir) / "manifest.json"
        md5s = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        programs = build_marker_programs(args.data_dir, level=args.level, top_n=args.top_n)
        out = export_programs(programs, args.out_dir, md5s, args.level, args.top_n)
        print(f"wrote {out} ({len(programs)} programs)")


if __name__ == "__main__":
    main()
