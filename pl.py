### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import scanpy as sc


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
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
def __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot, obsm_slot = None, uns_slot = None):
    adata_stored = adata.uns[uns_data_slot].copy()
    if obsm_slot is not None:
        adata_stored.obsm[obsm_slot] = adata.obsm[obsm_slot]
    if uns_slot is not None:
        adata_stored.obsm[uns_slot] = adata.obsm[uns_slot]
    adata = adata_stored
    return adata

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------
def pca(adata,
        *,
        plot_stored_gex_data=False,
        plot_stored_pax_data=False,
        **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.pca.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm['X_pca'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm['X_pca'].
    **kwargs
        Arguments to provide to the sc.pl.pca function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot='X_pca')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot='X_pca')
    sc.pl.pca(adata, **kwargs)

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
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm['X_umap'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm['X_umap'].
    **kwargs
        Arguments to provide to the sc.pl.umap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot='X_umap')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot='X_umap')
    sc.pl.umap(adata, **kwargs)

def tsne(adata,
         *,
         plot_stored_gex_data=False,
         plot_stored_pax_data=False,
         **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.tsne.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm['X_tsne'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm['X_tsne'].
    **kwargs
        Arguments to provide to the sc.pl.tsne function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot='X_tsne')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot='X_tsne')
    sc.pl.tsne(adata, **kwargs)

def diffmap(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.diffmap.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm['X_diffmap'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm['X_diffmap'].
    **kwargs
        Arguments to provide to the sc.pl.diffmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot='X_diffmap')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot='X_diffmap')
    sc.pl.diffmap(adata, **kwargs)

def draw_graph(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.draw_graph.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm['X_draw_graph_fa'] or adata.obsm['X_draw_graph_fr'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm['X_draw_graph_fa'] or adata.obsm['X_draw_graph_fr'].
    **kwargs
        Arguments to provide to the sc.pl.draw_graph function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_stored_gex_data or plot_stored_pax_data:
        if plot_stored_gex_data:
            uns_data_slot = 'gex_data'
        elif plot_stored_pax_data:
            uns_data_slot = 'pax_data'
        adata_stored = adata.uns[uns_data_slot].copy()

        if "X_draw_graph_fa" in adata_stored.obsm.keys():
            adata_stored.obsm["X_draw_graph_fa"] = adata.obsm["X_draw_graph_fa"]
        if "X_draw_graph_fr" in adata_stored.obsm.keys():
            adata_stored.obsm["X_draw_graph_fr"] = adata.obsm["X_draw_graph_fr"]
        adata = adata_stored

    sc.pl.draw_graph(adata, **kwargs)

def spatial(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.spatial.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.uns['spatial'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.uns['spatial'].
    **kwargs
        Arguments to provide to the sc.pl.spatial function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=None,
                                                       uns_slot='spatial')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot=None,
                                                       uns_slot='spatial')
    sc.pl.spatial(adata, **kwargs)

def embedding(adata,
              *,
              basis,
              plot_stored_gex_data=False,
              plot_stored_pax_data=False,
              **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.embedding.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    basis
        The name of the represenation in adata.obsm that should be used for plotting.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm[basis].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm[basis].
    **kwargs
        Arguments to provide to the sc.pl.embedding function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=basis)
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot=basis)
    sc.pl.embedding(adata, basis, **kwargs)

def embedding_density(adata,
                      *,
                      basis='umap',
                      plot_stored_gex_data=False,
                      plot_stored_pax_data=False,
                      **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.embedding_density.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    basis (default: 'umap')
        The name of the represenation in adata.obsm that should be used for plotting.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'] on adata.obsm[basis].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'] on adata.obsm[basis].
    **kwargs
        Arguments to provide to the sc.pl.embedding_density function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=basis)
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='pax_data',
                                                       obsm_slot=basis)
    sc.pl.embedding_density(adata, basis, **kwargs)

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
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.scatter function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
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
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.heatmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.heatmap(adata,**kwargs)

def dotplot(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.dotplot.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.dotplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.dotplot(adata,**kwargs)

def tracksplot(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.tracksplot.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.tracksplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.tracksplot(adata,**kwargs)

def violin(adata,
           *,
           plot_stored_gex_data=False,
           plot_stored_pax_data=False,
           **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.violin.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.violin function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.violin(adata,**kwargs)

def stacked_violin(adata,
                   *,
                   plot_stored_gex_data=False,
                   plot_stored_pax_data=False,
                   **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.stacked_violin.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.stacked_violin function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.stacked_violin(adata,**kwargs)

def matrixplot(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.matrixplot.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.matrixplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.matrixplot(adata,**kwargs)

def clustermap(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.clustermap.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.clustermap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.clustermap(adata,**kwargs)

def ranking(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.ranking.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.ranking function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.ranking(adata,**kwargs)

def dendrogram(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.dendrogram.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_stored_gex_data (default: False)
        Plot adata.uns['gex_data'].
    plot_stored_pax_data (default: False)
        Plot adata.uns['pax_data'].
    **kwargs
        Arguments to provide to the sc.pl.dendrogram function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
    elif(plot_stored_pax_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='pax_data')
    sc.pl.dendrogram(adata,**kwargs)
