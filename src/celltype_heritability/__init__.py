"""Utilities for mapping GWAS heritability to brain cell types (Paper 1).

Data-free implementations of the S-LDSC / MAGMA / sc-linker plumbing used by
this project. Production runs use the official tools; these modules exist for
testing, prototyping, and transparent review.
"""

from . import annotations, ldscore, magma, sclinker, simulate, sldsc

__all__ = ["annotations", "ldscore", "magma", "sclinker", "simulate", "sldsc"]
__version__ = "0.1.0"
