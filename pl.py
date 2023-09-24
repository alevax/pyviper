### ---------- IMPORT DEPENDENCIES ----------
# import numpy as np
# from scipy.stats import ttest_1samp
# import os
# import shutil
import pandas as pd
# import anndata
import scanpy as sc
# from pyther_classes import *
# from pyther_narnea import *


### ---------- EXPORT LIST ----------
__all__ = []
# __all__ = ['umap',
#            'scatter',
#            'heatmap',
#            'dotplot',
#            'tracksplot',
#            'violin',
#            'stacked_violin',
#            'matrixplot',
#            'clustermap',
#            'ranking',
#            'dendrogram']

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------
def umap(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.umap.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.umap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.umap(adata, **kwargs)

def scatter(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.scatter.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.scatter function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.scatter(adata,**kwargs)

def heatmap(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.heatmap.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.heatmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.heatmap(adata,**kwargs)

def dotplot(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.dotplot(adata,**kwargs)

def tracksplot(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.tracksplot(adata,**kwargs)

def violin(adata,
           *,
           plot_stored_gex_data=False,
           plot_stored_pax_data=False,
           **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.violin(adata,**kwargs)

def stacked_violin(adata,
                   *,
                   plot_stored_gex_data=False,
                   plot_stored_pax_data=False,
                   **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.stacked_violin(adata,**kwargs)

def matrixplot(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.matrixplot(adata,**kwargs)

def clustermap(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.clustermap(adata,**kwargs)

def ranking(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.ranking(adata,**kwargs)

def dendrogram(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.dendrogram(adata,**kwargs)
