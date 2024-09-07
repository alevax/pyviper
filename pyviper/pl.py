### ---------- IMPORT DEPENDENCIES ----------
import scanpy as sc
import pandas as pd
import seaborn
from ._pl_sns_heatmap import _completemap
from ._pl_sns_heatmap import _mrs

### ---------- EXPORT LIST ----------
__all__ = []

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
def __get_pax_params(kwargs):
    pax_kwargs = kwargs.copy()
    if 'cmap' not in kwargs:
        pax_kwargs['cmap'] = 'RdBu_r'
    if 'vcenter' in kwargs:
        pax_kwargs['vcenter'] = 0
    return pax_kwargs

def __get_gex_param(kwargs):
    if 'cmap' in kwargs:
        cmap = kwargs['cmap']
    else:
        cmap = 'viridis'

def __parse_color(adata, kwargs):
    adata_vars = list(adata.var_names.values) + list(adata.obs.columns.values)
    kwargs = kwargs.copy()
    if 'color' in kwargs:
        color = kwargs['color']
        color = [c for c in color if c in adata_vars]
        kwargs['color'] = color
    return kwargs

def __parse_var_names(adata, kwargs):
    adata_vars = list(adata.var_names.values)
    kwargs = kwargs.copy()
    if 'var_names' in kwargs:
        var_names = kwargs['var_names']
        var_names = [v for v in var_names if v in adata_vars]
        kwargs['var_names'] = var_names
    return kwargs

def __parse_keys(adata, kwargs):
    adata_vars = list(adata.var_names.values) + list(adata.obs.columns.values)
    kwargs = kwargs.copy()
    if 'keys' in kwargs:
        keys = kwargs['keys']
        keys = [k for k in keys if k in adata_vars]
        kwargs['keys'] = keys
    return kwargs

def __get_stored_uns_data_and_prep_to_plot(
    adata,
    uns_data_slot,
    obsm_slot = None,
    uns_slot = None
):
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
        plot_pax_data=False,
        plot_gex_data=False,
        **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.pca.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot VIPER stored in adata on adata.obsm['X_pca'].
    plot_gex_data : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_pca'].
    **kwargs
        Arguments to provide to the sc.pl.pca function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data is True:
        pax_kwargs = __get_pax_params(kwargs)
        pax_kwargs = __parse_color(adata, pax_kwargs)
        sc.pl.pca(adata, **pax_kwargs)
    if plot_gex_data is True:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot='X_pca'
        )
        kwargs = __parse_color(adata, kwargs)
        sc.pl.pca(adata, **kwargs)

def umap(adata,
         *,
         plot_pax_data=True,
         plot_gex_data=False,
         **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.umap.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot VIPER stored in adata on adata.obsm['X_umap'].
    plot_gex_data : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_umap'].
    **kwargs
        Arguments to provide to the sc.pl.pca function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data is True:
        pax_kwargs = __get_pax_params(kwargs)
        pax_kwargs = __parse_color(adata, pax_kwargs)
        sc.pl.umap(adata, **pax_kwargs)
    if plot_gex_data is True:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot='X_umap'
        )
        kwargs = __parse_color(adata, kwargs)
        sc.pl.umap(adata, **kwargs)

def tsne(adata,
         *,
         plot_pax_data=True,
         plot_gex_data=False,
         **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.tsne.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot VIPER stored in adata on adata.obsm['X_tsne'].
    plot_gex_data : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_tsne'].
    **kwargs
        Arguments to provide to the sc.pl.tsne function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data is True:
        pax_kwargs = __get_pax_params(kwargs)
        pax_kwargs = __parse_color(adata, pax_kwargs)
        sc.pl.tsne(adata, **pax_kwargs)
    if plot_gex_data is True:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot='X_tsne'
        )
        kwargs = __parse_color(adata, kwargs)
        sc.pl.tsne(adata, **kwargs)

def diffmap(adata,
            *,
            plot_pax_data=True,
            plot_gex_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.diffmap.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot VIPER stored in adata on adata.obsm['X_diffmap'].
    plot_gex_data : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_diffmap'].
    **kwargs
        Arguments to provide to the sc.pl.diffmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data is True:
        pax_kwargs = __get_pax_params(kwargs)
        pax_kwargs = __parse_color(adata, pax_kwargs)
        sc.pl.diffmap(adata, **pax_kwargs)
    if plot_gex_data is True:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot='X_diffmap'
        )
        kwargs = __parse_color(adata, kwargs)
        sc.pl.diffmap(adata, **kwargs)

def draw_graph(adata,
               *,
               plot_pax_data=True,
               plot_gex_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.draw_graph.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot VIPER stored in adata on adata.obsm['X_draw_graph_fa'] or adata.obsm['X_draw_graph_fr'].
    plot_gex_data : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_draw_graph_fa'] or adata.obsm['X_draw_graph_fr'].
    **kwargs
        Arguments to provide to the sc.pl.draw_graph function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data is True:
        pax_kwargs = __get_pax_params(kwargs)
        pax_kwargs = __parse_color(adata, pax_kwargs)
        sc.pl.draw_graph(adata, **pax_kwargs)

    if plot_gex_data is True:
        uns_data_slot = 'gex_data'
        adata_stored = adata.uns[uns_data_slot].copy()

        if "X_draw_graph_fa" in adata_stored.obsm.keys():
            adata_stored.obsm["X_draw_graph_fa"] = adata.obsm["X_draw_graph_fa"]
        if "X_draw_graph_fr" in adata_stored.obsm.keys():
            adata_stored.obsm["X_draw_graph_fr"] = adata.obsm["X_draw_graph_fr"]
        adata = adata_stored

        kwargs = __parse_color(adata, kwargs)
        sc.pl.draw_graph(adata, **kwargs)

def spatial(adata,
            *,
            plot_pax_data=True,
            plot_gex_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.spatial.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata on adata.uns['spatial'].
    plot_gex_data : default: False
        Plot adata.uns['gex_data'] on adata.uns['spatial'].
    **kwargs
        Arguments to provide to the sc.pl.spatial function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_color(adata, kwargs)
        sc.pl.spatial(adata, **pax_kwargs)
    if(plot_gex_data is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=None,
                                                       uns_slot='spatial')
        kwargs = __parse_color(adata, kwargs)
        sc.pl.spatial(adata, **kwargs)


def embedding(adata,
              *,
              basis,
              plot_pax_data=True,
              plot_gex_data=False,
              **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.embedding.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    basis
        The name of the represenation in adata.obsm that should be used for plotting.
    plot_pax_data : default: True
        Plot adata on adata.obsm[basis].
    plot_gex_data : default: True
        Plot adata.uns['gex_data'] on adata.obsm[basis].
    **kwargs
        Arguments to provide to the sc.pl.embedding function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_color(adata, kwargs)
        sc.pl.embedding(adata, basis, **pax_kwargs)
    if plot_gex_data:
        kwargs = __parse_color(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot=basis
        )
        sc.pl.embedding(adata, basis, **kwargs)

def embedding_density(adata,
                      *,
                      basis='umap',
                      plot_pax_data=True,
                      plot_gex_data=False,
                      **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.embedding_density.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    basis : default: 'umap'
        The name of the represenation in adata.obsm that should be used for plotting.
    plot_pax_data : default: True
        Plot adata on adata.obsm[basis].
    plot_gex_data : default: False
        Plot adata.uns['gex_data'] on adata.obsm[basis].
    **kwargs
        Arguments to provide to the sc.pl.embedding_density function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        sc.pl.embedding_density(adata, basis, **kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=basis)
        sc.pl.embedding_density(adata, basis, **kwargs)

def scatter(adata,
            *,
            plot_pax_data=True,
            plot_gex_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.scatter.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.scatter function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_color(adata, kwargs)
        sc.pl.scatter(adata,**pax_kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_color(adata, kwargs)
        sc.pl.scatter(adata,**kwargs)

def heatmap(adata,
            *,
            plot_pax_data=True,
            plot_gex_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.heatmap.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.heatmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_var_names(adata, kwargs)
        sc.pl.heatmap(adata,**pax_kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.heatmap(adata,**kwargs)

def dotplot(adata,
            *,
            plot_pax_data=True,
            plot_gex_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.dotplot.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot VIPER stored in adata.
    plot_gex_data : default: False
        Plot gExpr stored in adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.dotplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data is True:
        pax_kwargs = kwargs.copy()
        if 'cmap' not in kwargs:
            pax_kwargs['cmap'] = 'Reds'
        pax_kwargs = __parse_var_names(adata, kwargs)
        sc.pl.dotplot(adata, **pax_kwargs)
    if plot_gex_data is True:
        gex_kwargs = kwargs.copy()
        if 'cmap' not in kwargs:
            gex_kwargs['cmap'] = 'Greens'
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.dotplot(adata, **gex_kwargs)

def tracksplot(adata,
               *,
               plot_pax_data=True,
               plot_gex_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.tracksplot.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.tracksplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_var_names(adata, kwargs)
        sc.pl.tracksplot(adata,**pax_kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.tracksplot(adata,**kwargs)

def violin(adata,
           *,
           plot_pax_data=True,
           plot_gex_data=False,
           **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.violin.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.violin function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_keys(adata, kwargs)
        sc.pl.violin(adata,**pax_kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_keys(adata, kwargs)
        sc.pl.violin(adata,**kwargs)

def stacked_violin(adata,
                   *,
                   plot_pax_data=True,
                   plot_gex_data=False,
                   **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.stacked_violin.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.stacked_violin function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_var_names(adata, pax_kwargs)
        sc.pl.stacked_violin(adata,**pax_kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.stacked_violin(adata,**kwargs)

def matrixplot(adata,
               *,
               plot_pax_data=True,
               plot_gex_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.matrixplot.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.matrixplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        pax_kwargs = __parse_var_names(adata, pax_kwargs)
        sc.pl.matrixplot(adata,**pax_kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.matrixplot(adata,**kwargs)

def clustermap(adata,
               *,
               plot_pax_data=True,
               plot_gex_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.clustermap.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.clustermap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        sc.pl.clustermap(adata,**kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        sc.pl.clustermap(adata,**kwargs)

def ranking(adata,
            *,
            plot_pax_data=True,
            plot_gex_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.ranking.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.ranking function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        sc.pl.ranking(adata,**kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
        sc.pl.ranking(adata,**kwargs)

def dendrogram(adata,
               *,
               plot_pax_data=True,
               plot_gex_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.dendrogram.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax_data : default: True
        Plot adata.
    plot_gex_data : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.dendrogram function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax_data:
        sc.pl.dendrogram(adata,**kwargs)
    if plot_gex_data:
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')

def completemap(
    adata,
    var_names,
    cluster_column=None,
    show_gex_heatmap=False,
    plot_gex_norm=False,
    obs_metadata=None,
    gex_metadata=None,
    pax_metadata=None,
    h_clust_rows=True,
    h_clust_cols=False
):
    """\
    A function to display VIPER and gExpr along with multiple rows of metadata with samples organized by silhouette score within each cluster.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    var_names
        A list of genes/proteins to be visualized in the heatmap.
    cluster_column : default: None
        A column in adata.obs to sort the samples based on. Samples will be
        arranged from highest silhouette score to lowest silhouette score within
        each cluster.
    show_gex_heatmap : default: False
        Whether to plot var_names on an additional gExpr heatmap.
    plot_gex_norm : default: False
        Whether gExpr in the gex_heatmap and gex_metadata is normalized gExpr
        instead of scaled. Requires adata.uns['gex_data'].raw to contain
        unnormalized counts.
    obs_metadata : default: None
        A str or list of columns in adata.obs to visualize as an annotation
        above the heatmap.
    gex_metadata : default: None
        A str or list of vars in adata.uns['gex_data'] to visualize as an annotation above the heatmap.
    pax_metadata : default: None
        A str or list of vars in adata to visualize as an annotation above the heatmap.
    h_clust_rows : default: True
        Whether to hierarchically cluster the rows.
    h_clust_cols : default: False
        Whether to hierarchically cluster the columns. With False and cluster_column supplied, samples will be ordered by silhouette score.

    Returns
    -------
    A plot of :class:`ClusterGrid`.

    """
    _completemap(
        adata,
        var_names,
        cluster_column,
        show_gex_heatmap,
        plot_gex_norm,
        obs_metadata,
        gex_metadata,
        pax_metadata,
        h_clust_rows,
        h_clust_cols
    )

def mrs(
    adata,
    cluster_column,
    show_gex_heatmap=False,
    plot_gex_norm=False,
    obs_metadata=None,
    gex_metadata=None,
    pax_metadata=None,
    h_clust_rows=True,
    h_clust_cols=False,
    n_top_mrs=10,
    method='stouffer',
    top_mrs_list=None,
    mr_col=None
):
    """\
    A wrapper around pyviper.pl.completemap to visualize the top master regulators (MRs) from each cluster using the VIPER anndata object.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    var_names
        A list of genes/proteins to be visualized in the heatmap.
    cluster_column : default: None
        A column in adata.obs to sort the samples based on. Samples will be
        arranged from highest silhouette score to lowest silhouette score within
        each cluster.
    show_gex_heatmap : default: False
        Whether to plot var_names on an additional gExpr heatmap.
    plot_gex_norm : default: False
        Whether gExpr in the gex_heatmap and gex_metadata is normalized gExpr
        instead of scaled. Requires adata.uns['gex_data'].raw to contain
        unnormalized counts.
    obs_metadata : default: None
        A str or list of columns in adata.obs to visualize as an annotation
        above the heatmap.
    gex_metadata : default: None
        A str or list of vars in adata.uns['gex_data'] to visualize as an annotation above the heatmap.
    pax_metadata : default: None
        A str or list of vars in adata to visualize as an annotation above the heatmap.
    h_clust_rows : default: True
        Whether to hierarchically cluster the rows.
    h_clust_cols : default: False
        Whether to hierarchically cluster the columns. With False and cluster_column supplied, samples will be ordered by silhouette score.
    n_top_mrs : default: 10
        If top_mrs_list is None and mr_col is None, how many MRs per cluster to identify.
    method : default: 'stouffer'
        If top_mrs_list is None and mr_col is None, what method should be used to calculate the top MRs. Options include 'stouffer' (Stouffer signature), "mwu" (Mann-Whitney U-Test), 'spearman' (correlation of proteins with proximity to each cluster's center).
    top_mrs_list : default: None
        Whether to manually supply a list of MRs to be visualized.
    mr_col : default: None
        Whether to select a column from adata.var where MRs are marked, e.g. by running pyviper.tl.find_top_mrs.
    """

    _mrs(
        adata,
        cluster_column,
        show_gex_heatmap,
        plot_gex_norm,
        obs_metadata,
        gex_metadata,
        pax_metadata,
        h_clust_rows,
        h_clust_cols,
        n_top_mrs,
        method,
        top_mrs_list,
        mr_col
    )
