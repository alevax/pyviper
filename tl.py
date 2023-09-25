### ---------- IMPORT DEPENDENCIES ----------
# import numpy as np
# from scipy.stats import ttest_1samp
# import os
# import shutil
import pandas as pd
# import anndata
import scanpy as sc
from ._filtering_funcs import *
from ._filtering_funcs import _get_anndata_filtered_by_feature_group
# from pyther_classes import *
# from pyther_narnea import *
# import pathlib

### ---------- EXPORT LIST ----------
__all__ = []

# ------------------------ SCANPY TOOLS PYTHER WRAPPERS -----------------------
def pca(adata,
           *,
           layer = None,
           filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
           **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)

    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]
    return(adata)

def tsne(adata,
            *,
            layer = None,
            filter_by_feature_groups=None,  # ["tfs", "cotfs", "sig", "surf"],
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.tsne(adata_filt, **kwargs)
    adata.obsm["X_tsne"] = adata_filt.obsm["X_tsne"]
    return(adata)

def umap(adata,
            *,
            layer = None,
            filter_by_feature_groups=None,  # ["tfs", "cotfs", "sig", "surf"],
            svd_solver='arpack',
            n_neighbors=10,
            n_pcs=40,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    # Create a second anndata object
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # Move sc.pp.neighbors results to adata_filt
    # adata_filt.obsp["distances"] = adata.obsp["distances"]
    # adata_filt.obsp["connectivities"] = adata.obsp["connectivities"]
    # adata_filt.uns["neighbors"] = adata.uns["neighbors"]
    # Compute the UMAP
    sc.tl.umap(adata_filt, **kwargs)
    # Give the UMAP from the second anndata object to the original
    adata.obsm["X_umap"] = adata_filt.obsm["X_umap"]
    return(adata)

def draw_graph(adata,
                  *,
                  layer = None,
                  filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
                  **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.draw_graph(adata_filt, **kwargs)
    if "X_draw_graph_fa" in list(adata.obsm.keys()):
        adata.obsm["X_draw_graph_fa"] = adata_filt.obsm["X_draw_graph_fa"]
    if "X_draw_graph_fr" in list(adata.obsm.keys()):
        adata.obsm["X_draw_graph_fr"] = adata_filt.obsm["X_draw_graph_fr"]
    return(adata)

def diffmap(adata,
               *,
               layer = None,
               filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
               **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.diffmap(adata_filt, **kwargs)
    adata.obsm["X_diffmap"] = adata_filt.obsm["X_diffmap"]
    adata.uns['diffmap_evals'] = adata_filt.uns['diffmap_evals']
    return(adata)
