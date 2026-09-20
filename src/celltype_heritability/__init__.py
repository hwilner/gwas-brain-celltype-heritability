"""Utilities for mapping GWAS heritability to brain cell types (Paper 1).

Data-free implementations of the S-LDSC / MAGMA / sc-linker plumbing used by
this project, plus loaders for the public ABC Atlas marker tables and GWAS
summary statistics. Production runs use the official tools; these modules
exist for testing, prototyping, and transparent review.
"""

from . import annotations, atlas, ldscore, magma, sclinker, simulate, sldsc, sumstats

__all__ = [
    "annotations",
    "atlas",
    "ldscore",
    "magma",
    "sclinker",
    "simulate",
    "sldsc",
    "sumstats",
]
__version__ = "0.1.0"
