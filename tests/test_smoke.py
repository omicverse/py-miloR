"""Smoke tests that exercise the import surface and basic construction.

The milo pipeline needs a real MuData/AnnData with k-NN + sample metadata,
which is heavy to fabricate here — integration testing lives upstream in
omicverse. These tests only verify the standalone package imports and
constructs the core objects cleanly.
"""
import numpy as np
import anndata as ad
import pandas as pd
import pytest


def test_import_surface():
    import milor_py
    assert hasattr(milor_py, "Milo")
    assert hasattr(milor_py, "calcNormFactors")
    assert milor_py.__version__


def test_milo_constructs():
    from milor_py import Milo
    m = Milo()
    assert m is not None


def test_calc_norm_factors_runs_on_simple_counts():
    """TMM normalisation on a tiny synthetic count matrix."""
    from milor_py import calcNormFactors
    rng = np.random.default_rng(0)
    # 200 genes × 5 samples of Poisson counts
    counts = rng.poisson(50.0, size=(200, 5)).astype(np.float64)
    nf = calcNormFactors(counts, method="TMM")
    nf = np.asarray(nf)
    assert nf.shape == (5,)
    assert np.all(np.isfinite(nf))
    # All factors should be reasonably close to 1.0 for iid Poisson columns
    assert np.all(nf > 0.1) and np.all(nf < 10.0)
