### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import scanpy as sc
from ._filtering_funcs import *
from ._filtering_funcs import _get_anndata_filtered_by_feature_group

### ---------- EXPORT LIST ----------
__all__ = []

# ------------------------ SCANPY TOOLS PYTHER WRAPPERS -----------------------
def pca(adata,
        *,
        layer=None,
        filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
        **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)

    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]

def tsne(adata,
            *,
            layer = None,
            filter_by_feature_groups=["tfs", "cotfs"],  # ["tfs", "cotfs", "sig", "surf"],
            **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.tsne(adata_filt, **kwargs)
    adata.obsm["X_tsne"] = adata_filt.obsm["X_tsne"]

def umap(adata,
         *,
         layer=None,
         filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
         **kwargs):
    # Create a second anndata object
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    # Compute the UMAP
    sc.tl.umap(adata_filt, **kwargs)
    # Give the UMAP from the second anndata object to the original
    adata.obsm["X_umap"] = adata_filt.obsm["X_umap"]

def draw_graph(adata,
               *,
               layer=None,
               filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
               **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.draw_graph(adata_filt, **kwargs)
    if "X_draw_graph_fa" in list(adata_filt.obsm.keys()):
        adata.obsm["X_draw_graph_fa"] = adata_filt.obsm["X_draw_graph_fa"]
    if "X_draw_graph_fr" in list(adata_filt.obsm.keys()):
        adata.obsm["X_draw_graph_fr"] = adata_filt.obsm["X_draw_graph_fr"]

def diffmap(adata,
            *,
            layer=None,
            filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
            **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.diffmap(adata_filt, **kwargs)
    adata.obsm["X_diffmap"] = adata_filt.obsm["X_diffmap"]
    adata.uns['diffmap_evals'] = adata_filt.uns['diffmap_evals']

def leiden(adata,
           *,
           key_added='leiden',
           layer=None,
           filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
           **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.leiden(adata_filt, **kwargs, key_added = key_added)
    adata.obs[key_added] = adata_filt.obs[key_added]
    adata.uns['leiden'] = adata_filt.uns['leiden']

def louvain(adata,
            *,
            key_added='louvain',
            layer=None,
            filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
            **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.louvain(adata_filt, **kwargs)
    adata.obs[key_added] = adata_filt.obs[key_added]

def dendrogram(adata,
               *,
               groupby,
               key_added=None,
               layer=None,
               filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
               **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    if key_added is None:
        key_added = f'dendrogram_{groupby}'
    sc.tl.dendrogram(adata_filt, groupby, **kwargs, key_added = key_added)
    adata.uns[key_added] = adata_filt.uns[key_added]

def dpt(adata,
        *,
        layer=None,
        filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
        **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.dpt(adata_filt, **kwargs)
    adata.obs['dpt_pseudotime'] = adata_filt.obs['dpt_pseudotime']
    adata.obs['dpt_groups'] = adata_filt.obs['dpt_groups']

def paga(adata,
         groups="leiden",
         *,
         layer=None,
         filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
         **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.paga(adata_filt, groups, **kwargs)
    adata.uns['paga'] = adata_filt.uns['paga']
    adata.uns[str(groups) + '_sizes'] = adata_filt.uns[str(groups) + '_sizes']