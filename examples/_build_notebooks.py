"""Build and execute the Milo example notebooks.

Produces:
  - tutorial_omicverse.ipynb       — Haber 2017 DA via ov.single.DCT(method='milopy')
  - tutorial_standalone.ipynb      — Same analysis via from milor_py import Milo

Both use the Haber et al. 2017 (Nature) intestinal epithelium dataset
with Control vs Salmonella, following the omicverse t_deg_single
tutorial "DCT with milo" section.
"""
import os
import nbformat as nbf
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = "/scratch/users/steorra/analysis/26_omic_protocol/data/haber_2017_regions.h5ad"


# --------------------------------------------------------------------- #
# Tutorial via omicverse (the canonical entrypoint)
# --------------------------------------------------------------------- #
def _tut_omicverse():
    nb = nbf.v4.new_notebook()
    c = nb.cells

    c.append(nbf.v4.new_markdown_cell("""\
# Milo tutorial via **omicverse** (`ov.single.DCT(method='milopy')`)

Differential abundance testing on single-cell data using the Milo
k-NN-neighbourhood framework (Dann et al., *Nature Biotechnology* 2022).
This notebook is the "DCT with milo" section of the omicverse
`t_deg_single` tutorial, adapted as a standalone example.

**Data**: Haber et al. 2017 mouse intestinal epithelium (Control vs
Salmonella, ~5 k cells × 15 k genes).

**Backend**: this walks through the `ov.single.DCT` wrapper, which
dispatches to the canonical `omicverse.single.Milo` class — the same
code this repo mirrors as `from milor_py import Milo`.
"""))

    c.append(nbf.v4.new_code_cell(f"""\
import warnings; warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import omicverse as ov

ov.plot_set()
adata = ov.read({DATA_PATH!r})
print(adata)"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 1. Subset to Control vs Salmonella"""))

    c.append(nbf.v4.new_code_cell("""\
adata = adata[adata.obs['condition'].isin(['Control', 'Salmonella'])].copy()
print(adata.obs['condition'].value_counts())"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 2. Preprocess, batch-correct, build the k-NN graph

Milo's index sampling + neighbourhood expansion both operate on the
k-NN graph, so a clean batch-corrected PC space matters."""))

    c.append(nbf.v4.new_code_cell("""\
adata = ov.pp.preprocess(adata, mode='shiftlog|pearson',
                          n_HVGs=2000, target_sum=50*1e4)
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
ov.single.batch_correction(adata, batch_key='batch',
                            methods='harmony', n_pcs=50)
ov.pp.neighbors(adata, n_neighbors=15, n_pcs=50,
                use_rep='X_pca_harmony')
ov.pp.umap(adata)"""))

    c.append(nbf.v4.new_code_cell("""\
ov.pl.embedding(adata, basis='X_umap',
                color=['batch', 'cell_label'], frameon='small')"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 3. Run differential abundance with `ov.single.DCT`

`method='milopy'` routes through `omicverse.single.Milo` (the same
class this repo exports as `from milor_py import Milo`). Internally
it runs `make_nhoods → count_nhoods → da_nhoods → build_nhood_graph`
in one step."""))

    c.append(nbf.v4.new_code_cell("""\
dct = ov.single.DCT(
    adata, condition='condition',
    ctrl_group='Control', test_group='Salmonella',
    cell_type_key='cell_label', method='milopy',
    sample_key='batch', use_rep='X_pca_harmony',
)
dct.run()"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 4. Inspect DA results

`dct.mdata['milo']` is an AnnData keyed by *neighbourhood* (not cell).
Its `.var` carries `logFC`, `PValue`, `SpatialFDR` plus
`nhood_annotation` (majority-vote cell type) and
`nhood_annotation_frac` (fraction of cells carrying that label)."""))

    c.append(nbf.v4.new_code_cell("""\
milo_var = dct.mdata['milo'].var
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].hist(milo_var['PValue'], bins=50, color='#4c78a8')
axes[0].set_xlabel('P-Value'); axes[0].set_title('Null p-value uniformity check')
axes[1].scatter(milo_var['logFC'], -np.log10(milo_var['PValue']),
                 s=8, alpha=0.6)
axes[1].set_xlabel('logFC'); axes[1].set_ylabel('-log10(p)')
axes[1].set_title('Volcano')
plt.tight_layout(); plt.show()"""))

    c.append(nbf.v4.new_markdown_cell("""\
Many neighbourhoods span multiple cell types. We call those whose
dominant label covers < 60 % of cells "Mixed"."""))

    c.append(nbf.v4.new_code_cell("""\
res = dct.get_results(mix_threshold=0.6)
res.head()"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 5. Visualise neighbourhood DA"""))

    c.append(nbf.v4.new_code_cell("""\
color_dict = dict(zip(adata.obs['cell_label'].cat.categories,
                       ov.pl.green_color[:4] + ov.pl.purple_color))
color_dict['Mixed'] = '#c2c2c2'
dct.model.plot_da_beeswarm(dct.mdata, alpha=0.1, palette=color_dict)
plt.xticks(fontsize=10); plt.yticks(fontsize=10); plt.tight_layout()"""))

    c.append(nbf.v4.new_markdown_cell("""\
## Summary

- The `DCT(method='milopy')` wrapper runs the full Milo pipeline in
  one call, suitable for integration with the rest of the omicverse
  single-cell analysis stack.
- The object `dct.model` is an instance of the exact same Milo class
  this repo ships as `from milor_py import Milo`. See
  `tutorial_standalone.ipynb` for the direct-API version."""))

    out = os.path.join(HERE, "tutorial_omicverse.ipynb")
    with open(out, "w") as f:
        nbf.write(nb, f)
    return out


# --------------------------------------------------------------------- #
# Tutorial via the standalone package
# --------------------------------------------------------------------- #
def _tut_standalone():
    nb = nbf.v4.new_notebook()
    c = nb.cells

    c.append(nbf.v4.new_markdown_cell("""\
# Milo tutorial via **`milor_py`** standalone

Same analysis as `tutorial_omicverse.ipynb` but driving the Milo
class directly — no `omicverse` imports required. Useful if you only
want differential abundance testing without the full omicverse stack.

**Data**: Haber et al. 2017 mouse intestinal epithelium (Control vs
Salmonella).

The Milo pipeline here is the canonical four-step workflow:

1. `make_nhoods` — sample index cells and expand each into its
   k-NN neighbourhood
2. `count_nhoods` — build the (n_nhoods × n_samples) count matrix
3. `da_nhoods` — fit a QL-NB GLM per neighbourhood, return `logFC` +
   `PValue` + spatial-FDR
4. `build_nhood_graph` — overlap-based neighbourhood graph for the
   beeswarm / graph plots
"""))

    c.append(nbf.v4.new_code_cell(f"""\
import warnings; warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad

from milor_py import Milo

adata = ad.read_h5ad({DATA_PATH!r})
adata = adata[adata.obs['condition'].isin(['Control', 'Salmonella'])].copy()
print(adata)
print(adata.obs['condition'].value_counts())"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 1. Preprocess + k-NN graph

Standard scanpy pipeline (no omicverse): normalise, PCA, neighbours,
UMAP. The k-NN graph is what Milo operates on in step 2."""))

    c.append(nbf.v4.new_code_cell("""\
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3',
                             layer=None)
adata = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)"""))

    c.append(nbf.v4.new_code_cell("""\
sc.pl.umap(adata, color=['batch', 'cell_label'], frameon=False,
            wspace=0.3, show=False)
plt.show()"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 2. Milo pipeline — direct API

The four core calls. Inputs and outputs are written into a MuData
object so the per-neighbourhood results live alongside the per-cell
AnnData."""))

    c.append(nbf.v4.new_code_cell("""\
m = Milo()

# Wrap the AnnData in a MuData with slots ``rna`` (input) + ``milo``
# (results will be written into this slot).
mdata = m.load(adata)

# 1) Sample 10% of cells as neighbourhood index points. Neighbourhood
#    size is inherited from the ``sc.pp.neighbors(n_neighbors=...)``
#    call above.
m.make_nhoods(mdata, prop=0.1)

# 2) Build the (n_nhoods × n_samples) count matrix; sample labels come from `batch`.
mdata = m.count_nhoods(mdata, sample_col='batch')

# 3) Fit the QL-NB GLM per neighbourhood.
#    ``model_contrasts`` selects the Salmonella − Control contrast.
m.da_nhoods(
    mdata,
    design='~condition',
    model_contrasts='condition[Salmonella]-condition[Control]',
)

# 4) Overlap graph between neighbourhoods for the beeswarm plot.
m.build_nhood_graph(mdata, basis='X_umap')

print(mdata['milo'])
print(mdata['milo'].var[['logFC', 'PValue', 'SpatialFDR']].head())"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 3. Annotate neighbourhoods with majority-vote cell type"""))

    c.append(nbf.v4.new_code_cell("""\
m.annotate_nhoods(mdata, anno_col='cell_label')
# Mark low-confidence neighbourhoods (< 60 % of cells sharing the
# majority label) as ``Mixed`` — extend the categorical first.
anno = mdata['milo'].var['nhood_annotation']
if hasattr(anno, 'cat') and 'Mixed' not in anno.cat.categories:
    anno = anno.cat.add_categories(['Mixed'])
mix_mask = mdata['milo'].var['nhood_annotation_frac'] < 0.6
anno[mix_mask] = 'Mixed'
mdata['milo'].var['nhood_annotation'] = anno
mdata['milo'].var[['nhood_annotation',
                    'nhood_annotation_frac',
                    'logFC', 'SpatialFDR']].head()"""))

    c.append(nbf.v4.new_markdown_cell("""\
## 4. Diagnostics + visualisation"""))

    c.append(nbf.v4.new_code_cell("""\
milo_var = mdata['milo'].var
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].hist(milo_var['PValue'], bins=50, color='#4c78a8')
axes[0].set_xlabel('P-Value'); axes[0].set_title('Null p-value uniformity')
axes[1].scatter(milo_var['logFC'], -np.log10(milo_var['PValue']),
                 s=8, alpha=0.6)
axes[1].set_xlabel('logFC'); axes[1].set_ylabel('-log10(p)')
axes[1].set_title('Volcano')
plt.tight_layout(); plt.show()"""))

    c.append(nbf.v4.new_code_cell("""\
m.plot_da_beeswarm(mdata, alpha=0.1)
plt.xticks(fontsize=10); plt.yticks(fontsize=10); plt.tight_layout()"""))

    c.append(nbf.v4.new_markdown_cell("""\
## Summary

The standalone API mirrors the miloR R workflow 1-for-1:

```python
from milor_py import Milo
m = Milo()
mdata = m.load(adata)            # AnnData → MuData with 'rna' + 'milo' slots
m.make_nhoods(mdata, prop=0.1)
mdata = m.count_nhoods(mdata, sample_col='batch')
m.da_nhoods(mdata, design='~condition',
             model_contrasts='condition[Salmonella]-condition[Control]')
m.build_nhood_graph(mdata, basis='X_umap')
m.annotate_nhoods(mdata, anno_col='cell_label')
```

The class you just used is the **same code path** as
`ov.single.DCT(method='milopy')` in the omicverse tutorial — omicverse
is the upstream development home, and this repo mirrors the Milo class
with zero functional diff."""))

    out = os.path.join(HERE, "tutorial_standalone.ipynb")
    with open(out, "w") as f:
        nbf.write(nb, f)
    return out


def _execute(path):
    nb = nbf.read(path, as_version=4)
    client = NotebookClient(nb, timeout=1800, kernel_name="python3",
                             resources={"metadata": {"path": HERE}})
    client.execute()
    with open(path, "w") as f:
        nbf.write(nb, f)
    print(f"executed {path}")


if __name__ == "__main__":
    p1 = _tut_omicverse()
    p2 = _tut_standalone()
    _execute(p1)
    _execute(p2)
