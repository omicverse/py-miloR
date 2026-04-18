# milor-py

A **pure-Python re-implementation of Milo** (Dann et al., *Nature Biotech.* 2022) for single-cell differential abundance testing on k-nearest-neighbour neighbourhoods.

- AnnData / MuData-native — drop-in for the scanpy ecosystem
- **No `rpy2`**, no R install, no edgeR dependency — TMM normalisation is implemented directly in NumPy + SciPy
- Same API as the R [miloR](https://github.com/MarioniLab/miloR) workflow (`make_nhoods` → `count_nhoods` → `da_nhoods` → `build_nhood_graph`)

> This is a **standalone mirror** of the canonical implementation that lives in [`omicverse`](https://github.com/Starlitnightly/omicverse) (`omicverse.single.Milo` / `omicverse/single/_milo_dev.py`). All algorithmic work is developed upstream in omicverse and synced here for users who want Milo without the full omicverse stack.

## Install

```bash
pip install milor-py
```

## Quick-start

```python
import anndata as ad
import scanpy as sc
from milor_py import Milo

adata = ad.read_h5ad("mydata.h5ad")       # cells × genes
sc.pp.neighbors(adata, n_neighbors=30)    # build the kNN graph first

m = Milo()

# 1) Sample 10% of cells as neighbourhood "index" cells and expand each
#    index into a k-NN neighbourhood.
adata = m.make_nhoods(adata, prop=0.1, k=30)

# 2) Count how many cells of each sample are in each neighbourhood,
#    producing an (n_nhoods × n_samples) count matrix in
#    ``adata.uns['nhood_adata']``.
adata = m.count_nhoods(adata, sample_col="sample")

# 3) Fit the quasi-likelihood negative-binomial GLM per neighbourhood
#    (implemented with statsmodels QL-F; matches miloR's edgeR step).
m.da_nhoods(adata, design="~condition")

# 4) Build the neighbourhood graph for visualisation.
m.build_nhood_graph(adata)
```

Results are written back into the AnnData / MuData object:

| Slot | Contents |
|---|---|
| `adata.obsm['nhoods']` | sparse cell × neighbourhood membership matrix |
| `adata.uns['nhood_adata']` | AnnData of (n_nhoods × n_samples) counts + per-neighbourhood DA results |
| `adata.uns['nhood_adata'].var['logFC']` / `['PValue']` / `['SpatialFDR']` | per-neighbourhood log-fold-change, raw p-value, spatial-FDR-corrected p-value |
| `adata.uns['nhood_adata'].obsm['X_nhood_graph']` | neighbourhood-graph coordinates for `plot_nhood_graph` |

## What's included

The `Milo` class mirrors the miloR Python bindings:

| Method | Purpose |
|---|---|
| `make_nhoods` | sample index cells and expand into k-NN neighbourhoods |
| `count_nhoods` | build the neighbourhood × sample count matrix |
| `da_nhoods` | fit QL-NB GLM per neighbourhood, return DA results |
| `annotate_nhoods` | assign each neighbourhood a majority-vote label |
| `annotate_nhoods_continuous` | continuous feature summary per neighbourhood |
| `add_covariate_to_nhoods_var` | pull a per-cell covariate into `nhood_adata.var` |
| `build_nhood_graph` | overlap-based neighbourhood graph for visualisation |
| `add_nhood_expression` | mean expression matrix per neighbourhood |
| `plot_nhood_graph` | coloured-by-logFC neighbourhood graph |
| `plot_nhood` / `plot_da_beeswarm` / `plot_nhood_counts_by_cond` | standard miloR-style diagnostics |

TMM normalisation (`calcNormFactors`) is also exposed at the top level — it's a direct NumPy port of edgeR's algorithm, useful independently of Milo.

## Notebooks

Two executed tutorials live under [`examples/`](examples/) — both run
the same Haber et al. 2017 (*Nature*) mouse-intestine dataset with
Control vs Salmonella and produce the same DA results.

| Notebook | Backend |
|---|---|
| [`examples/tutorial_omicverse.ipynb`](examples/tutorial_omicverse.ipynb) | `ov.single.DCT(method='milopy')` — the canonical entrypoint when you already use omicverse. |
| [`examples/tutorial_standalone.ipynb`](examples/tutorial_standalone.ipynb) | `from milor_py import Milo` — direct `make_nhoods → count_nhoods → da_nhoods → build_nhood_graph` pipeline, no omicverse required. |

Either notebook drives the identical `Milo` class; omicverse is the
upstream development home and this repo mirrors it.

## Relationship to omicverse

Developed **upstream** in [`omicverse`](https://github.com/Starlitnightly/omicverse):

- Canonical implementation: `omicverse.single.Milo` (`omicverse/single/_milo_dev.py`)
- Standalone mirror (this repo): same code, same API, minus the omicverse packaging

If you already use omicverse, there is no reason to install this package separately — `omicverse.single.Milo` exposes the same class. This repo exists for users who want differential abundance testing without the full omicverse stack.

## Relationship to `milopy`

A separate Python port of Milo exists as [`milopy`](https://github.com/emdann/milopy) by the original Milo author. The two projects share the same scientific algorithm but are independent implementations with different API surfaces. `milor-py` stays closer to the miloR R API and is developed alongside the broader omicverse single-cell stack.

## Citation

If you use this package, please cite the original Milo paper:

> Dann, E., Henderson, N.C., Teichmann, S.A., Morgan, M.D. & Marioni, J.C. **Differential abundance testing on single-cell data using k-nearest neighbor graphs.** *Nature Biotechnology* 40, 245–253 (2022).

and acknowledge omicverse / this repo for the Python port.

## License

GNU GPLv3 — matches [`omicverse`](https://github.com/Starlitnightly/omicverse) upstream.
