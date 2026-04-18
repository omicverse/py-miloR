"""
milor_py: pure-Python Milo for AnnData-native differential abundance testing.

A standalone mirror of the canonical implementation that lives in
``omicverse.single.Milo`` / ``omicverse/single/_milo_dev.py``. This
repo exists for users who want Milo without pulling in the full
omicverse stack — all algorithmic work is developed upstream in
omicverse and synced here.

Quick-start
-----------
>>> from milor_py import Milo
>>> m = Milo()
>>> adata = m.make_nhoods(adata, prop=0.1, k=30)
>>> adata = m.count_nhoods(adata, sample_col="sample")
>>> m.da_nhoods(adata, design="~condition")
>>> m.build_nhood_graph(adata)
"""
from __future__ import annotations

from .milo import (
    Milo,
    calcNormFactors,
)

__version__ = "0.1.0"

__all__ = [
    "Milo",
    "calcNormFactors",
]
