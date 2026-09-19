import numpy as np
import pandas as pd

from celltype_heritability.annotations import (
    add_flank,
    binarize_annotation,
    build_celltype_annotations,
    table_to_bed,
)


def _snps():
    return pd.DataFrame({"chrom": ["1"] * 5, "pos": [0, 50, 100, 150, 200]})


def test_table_to_bed_normalizes_and_sorts():
    t = pd.DataFrame(
        {"chrom": ["1", "1"], "start": [100, 0], "end": [120, 40], "gene": ["b", "a"]}
    )
    bed = table_to_bed(t, name_col="gene")
    assert list(bed.columns) == ["chrom", "start", "end", "name", "score"]
    assert list(bed["start"]) == [0, 100]


def test_add_flank_clips_at_zero():
    bed = pd.DataFrame({"chrom": ["1"], "start": [10], "end": [20], "name": ["."], "score": [0.0]})
    out = add_flank(bed, 50)
    assert out.loc[0, "start"] == 0 and out.loc[0, "end"] == 70


def test_binarize_annotation_overlaps():
    bed = pd.DataFrame({"chrom": ["1"], "start": [40], "end": [110], "name": ["."], "score": [0.0]})
    flags = binarize_annotation(bed, _snps())
    assert flags.tolist() == [0, 1, 1, 0, 0]


def test_build_celltype_annotations_with_window():
    beds = {
        "ct_a": pd.DataFrame({"chrom": ["1"], "start": [45], "end": [55], "name": ["."], "score": [0.0]}),
        "ct_b": pd.DataFrame({"chrom": ["1"], "start": [180], "end": [190], "name": ["."], "score": [0.0]}),
    }
    annot = build_celltype_annotations(beds, _snps(), window=50)
    assert annot["ct_a"].tolist() == [1, 1, 1, 0, 0]
    assert annot["ct_b"].tolist() == [0, 0, 0, 1, 1]
