import json

import pandas as pd
import pytest

from celltype_heritability import atlas


def _fake_atlas_dir(tmp_path):
    # minimal taxonomy: 2 clusters (1 neuronal, 1 non-neuronal), 2 superclusters
    terms = pd.DataFrame(
        [
            ("CS1", "Neuron A", "CCN202210140_CLUS"),
            ("CS2", "Microglia A", "CCN202210140_CLUS"),
            ("CS10", "Excitatory", "CCN202210140_SUPC"),
            ("CS11", "Microglia", "CCN202210140_SUPC"),
            ("CS20", "VGLUT1", "CCN202210140_NEUR"),
            ("CS21", None, "CCN202210140_NEUR"),
        ],
        columns=["label", "name", "cluster_annotation_term_set_label"],
    )
    terms.to_csv(tmp_path / "cluster_annotation_term.csv", index=False)
    pd.DataFrame({"cluster_alias": [0, 1], "number_of_cells": [10, 20],
                  "label": ["SUB1", "SUB2"]}).to_csv(tmp_path / "cluster.csv", index=False)
    mem = pd.DataFrame(
        [
            ("CS1", "CCN202210140_CLUS", 0),
            ("CS2", "CCN202210140_CLUS", 1),
            ("CS20", "CCN202210140_NEUR", 0),
            ("CS21", "CCN202210140_NEUR", 1),
        ],
        columns=["cluster_annotation_term_label",
                 "cluster_annotation_term_set_label", "cluster_alias"],
    )
    mem.to_csv(tmp_path / "cluster_to_cluster_annotation_membership.csv", index=False)
    markers = {
        "CCN202210140_CLUS/CS1": ["E1", "E2", "E3"],
        "CCN202210140_CLUS/CS2": ["E4", "E5", "E6"],
        "CCN202210140_SUPC/CS10": ["E1", "E2"],
        "CCN202210140_SUPC/CS11": ["E4", "E5"],
        "metadata": {"version": "test"},
    }
    (tmp_path / "query_markers.n10.20240221800.json").write_text(json.dumps(markers))
    genes = pd.DataFrame(
        {
            "gene_identifier": [f"E{i}" for i in range(1, 7)],
            "gene_symbol": [f"G{i}" for i in range(1, 7)],
            "biotype": ["protein_coding"] * 5 + ["lncRNA"],
        }
    )
    genes.to_csv(tmp_path / "gene.csv", index=False)
    return tmp_path


def test_cluster_programs_respect_top_n_and_protein_coding(tmp_path):
    d = _fake_atlas_dir(tmp_path)
    progs = atlas.build_marker_programs(d, level="cluster", top_n=2)
    assert progs == {"Neuron A": ["G1", "G2"], "Microglia A": ["G4", "G5"]}
    # deterministic
    assert atlas.build_marker_programs(d, level="cluster", top_n=2) == progs


def test_supercluster_programs_use_term_names(tmp_path):
    d = _fake_atlas_dir(tmp_path)
    progs = atlas.build_marker_programs(d, level="supercluster", top_n=5)
    assert progs == {"Excitatory": ["G1", "G2"], "Microglia": ["G4", "G5"]}


def test_class_programs_split_neuronal_nonneuronal(tmp_path):
    d = _fake_atlas_dir(tmp_path)
    progs = atlas.build_marker_programs(d, level="class", top_n=5)
    assert set(progs) == {"Neuronal", "Non-neuronal"}
    assert progs["Neuronal"] == ["G1", "G2", "G3"]
    assert progs["Non-neuronal"] == ["G4", "G5"]  # G6 (lncRNA) excluded


def test_export_programs_writes_json_and_manifest(tmp_path):
    out = tmp_path / "out"
    path = atlas.export_programs(
        {"ct": ["G1"]}, out, {"f": "abc"}, level="cluster", top_n=1
    )
    payload = json.loads(path.read_text())
    assert payload["programs"] == {"ct": ["G1"]}
    man = pd.read_csv(out / "programs_cluster_top1_manifest.csv")
    assert man.loc[0, "n_genes"] == 1
