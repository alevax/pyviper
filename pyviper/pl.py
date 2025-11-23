### ---------- IMPORT DEPENDENCIES ----------
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn
from ._pl_sns_heatmap import _ss_heatmap, _mrs
from ._pl_vis_net import _vis_net
from anndata import AnnData
import matplotlib.pyplot as plt

### ---------- EXPORT LIST ----------
__all__ = []

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
def __get_pax_params(kwargs, cmap_pax):
    pax_kwargs = kwargs.copy()
    if 'cmap' not in kwargs:
        pax_kwargs['cmap'] = cmap_pax
    if 'vcenter' in kwargs:
        pax_kwargs['vcenter'] = 0
    return pax_kwargs

def __remove_missing_items(adata, kwargs):
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
        if isinstance(var_names, list):
            var_names = [v for v in var_names if v in adata_vars]
        elif isinstance(var_names, dict):
            for key in var_names:
                if isinstance(var_names[key], str):
                    var_names[key] = [var_names[key]]
                var_names[key] = [v for v in var_names[key] if v in adata_vars]
        else:
            raise ValueError("var_names must be list or dict.")
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

def __get_indices_for_B_same_order_as_A(A, B):
    # Step 1: Create a mapping from values in B to their indices
    b_indices = {value: idx for idx, value in enumerate(B)}

    # Step 2: Use the mapping to create the index order for B that matches A
    reorder_indices = np.array([b_indices[value] for value in A])

    return reorder_indices

def __get_stored_uns_data_and_prep_to_plot(
    adata,
    uns_data_slot,
    obsm_slot = None,
    uns_slot = None
):
    adata_stored = adata.uns[uns_data_slot].copy()
    order_adata_as_adata_stored = __get_indices_for_B_same_order_as_A(
        adata_stored.obs_names.values, adata.obs_names.values
    )

    if obsm_slot is not None:
        adata_stored.obsm[obsm_slot] =\
            adata.obsm[obsm_slot][order_adata_as_adata_stored,:]
    if uns_slot is not None:
        adata_stored.obsm[uns_slot] =\
            adata.obsm[uns_slot][order_adata_as_adata_stored,:]
    adata = adata_stored
    return adata

def __get_adata_comb(adata, my_colors, basis = "X_pca"):
    obs_comb = pd.DataFrame(index=adata.obs_names)
    X_comb = pd.DataFrame(index=adata.obs_names)

    gex_order_as_pax = __get_indices_for_B_same_order_as_A(
                    adata.obs_names,
                    adata.uns['gex_data'].obs_names
    )

    missing_colors = []

    for c in my_colors:
        if c in adata.obs.columns:
            obs_comb[c] = adata.obs[c]
        elif c in adata.uns['gex_data'].obs:
            obs_comb[c] = adata.uns['gex_data'].obs[c][
                gex_order_as_pax
            ]
        else:
            if c in adata.var_names:
                X_comb[c] = adata.to_df()[c]
            if c in adata.uns["gex_data"].var_names:
                X_comb[c + "_gExpr"] = adata.uns["gex_data"].to_df()[c][gex_order_as_pax]
                my_colors = my_colors + [c + "_gExpr"]

            if c not in adata.var_names and c not in adata.uns["gex_data"].var_names:
                missing_colors+=[c]

    if X_comb.shape[1] > 0:
        adata_comb = AnnData(X_comb)
        adata_comb.obs = obs_comb
    else:
        adata_comb = AnnData(obs=obs_comb)
    adata_comb.obsm[basis] = adata.obsm[basis]

    missing_colors = set(missing_colors)
    my_colors = [item for item in my_colors if item not in missing_colors]

    return adata_comb, my_colors

def __change_gexpr_cmap(fig1, cmap_gex = 'viridis'):
    for i in range(len(fig1.get_children())):
        child = fig1.get_children()[i]

        # Check if the ax is an instance of Axes
        if isinstance(child, plt.Axes):
            # Extract the title from the Axes object
            title = child.get_title()
            # Modify the cmap of gExpr plots
            if "_gExpr" in title:
                path_collection = child.get_children()[0]
                path_collection.cmap = plt.get_cmap(cmap_gex)
                cbar = fig1.get_children()[i+1]
                cbar.get_children()[1].cmap = plt.get_cmap(cmap_gex)

def __change_obs_cmap(fig1, adata, cmap_obs = 'inferno'):
    for i in range(len(fig1.get_children())):
        child = fig1.get_children()[i]
        # Check if the ax is an instance of Axes
        if isinstance(child, plt.Axes):
            # Extract the title from the Axes object
            title = child.get_title()
            # Modify the cmap of gExpr plots
            if title in adata.obs.columns.values:
                is_numeric = pd.api.types.is_numeric_dtype(adata.obs[title])
                if is_numeric:
                    path_collection = child.get_children()[0]
                    path_collection.cmap = plt.get_cmap(cmap_obs)
                    cbar = fig1.get_children()[i+1]
                    cbar.get_children()[1].cmap = plt.get_cmap(cmap_obs)

def _plot(
    plot_func,
    basis,
    adata,
    plot_pax,
    plot_gex,
    cmap_pax,
    cmap_gex,
    cmap_obs,
    kwargs):

    if 'color' in kwargs and isinstance(kwargs['color'], str):
        kwargs['color'] = [kwargs['color']]

    if plot_pax is True and plot_gex is True:
        adata_combo, my_colors = __get_adata_comb(
            adata,
            kwargs['color'],
            basis = basis
        )
        kwargs['color'] = my_colors
        if 'cmap' not in kwargs:
            kwargs['cmap'] = cmap_pax
            kwargs['return_fig'] = True
            fig1 = plot_func(adata_combo, **kwargs)
            __change_gexpr_cmap(fig1, cmap_gex)
            __change_obs_cmap(fig1, adata, cmap_obs)
        else:
            plot_func(adata_combo, **kwargs)
    elif plot_pax is True:
        pax_kwargs = __get_pax_params(kwargs, cmap_pax)
        pax_kwargs = __remove_missing_items(adata, pax_kwargs)
        if 'cmap' not in kwargs:
            pax_kwargs['return_fig'] = True
            fig1 = plot_func(adata, **pax_kwargs)
            __change_obs_cmap(fig1, adata, cmap_obs)
        else:
            plot_func(adata, **pax_kwargs)
    elif plot_gex is True:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot=basis
        )
        kwargs = __remove_missing_items(adata, kwargs)
        if 'cmap' not in kwargs:
            kwargs['cmap'] = cmap_gex
            kwargs['return_fig'] = True
            fig1 = plot_func(adata, **kwargs)
            __change_obs_cmap(fig1, adata, cmap_obs)
        else:
            fig1 = plot_func(adata, **kwargs)

def _combo_dotplot(
    adata,
    cmap_pax,
    cmap_gex,
    spacing_factor,
    kwargs
):
    kwargs = kwargs.copy()
    groupby = kwargs['groupby']
    if 'show' in kwargs:
        plot_show = kwargs['show']
    else:
        plot_show = True

    kwargs['show'] = False

    # Get n_dots
    proteins = [x for x in kwargs['var_names'] if x in adata.var_names]
    genes = [x for x in kwargs['var_names'] \
             if x in adata.uns['gex_data'].var_names]
    del kwargs['var_names']

    # Create a new figure with the proportional size
    n_dots_max = np.max([len(proteins), len(genes)])
    n_categories = len(np.unique(adata.obs[groupby]))
    fig_new, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(3 + n_dots_max*0.5*spacing_factor,
                 3 + n_categories*1*spacing_factor)
    )

    # if groupby not in adata.uns['gex_data'].obs.columns:
    #     adata.uns['gex_data'].obs[groupby] = \
    #         adata.obs[groupby][__get_indices_for_B_same_order_as_A(
    #             adata.uns['gex_data'].obs_names.values,
    #             adata.obs_names.values
    #     )]
    __add_groupby_to_adata_gex_obs(adata, kwargs)

    pax_kwargs = kwargs.copy()
    if 'cmap' not in pax_kwargs:
        pax_kwargs['cmap'] = cmap_pax

    sc.pl.dotplot(
        adata,
        var_names = proteins,
        ax = ax1,
        **pax_kwargs,
    )

    gex_kwargs = kwargs.copy()
    if 'cmap' not in gex_kwargs:
        gex_kwargs['cmap'] = cmap_gex

    sc.pl.dotplot(
        adata.uns['gex_data'],
        var_names = genes,
        ax = ax2,
        **gex_kwargs
    )

    for ax in fig_new.axes:
        if ax.get_title() == 'Mean expression\nin group':
            selected_ax = ax
            break
    selected_ax.set_title('Mean activity\nin group')
    title_obj = selected_ax.title
    fontsize = title_obj.get_fontsize()
    selected_ax.set_title('Mean activity\nin group', fontsize=fontsize-3)

    if plot_show: plt.show()

def _violin_grid(
    adata,
    n_cols = 4,
    w_spacing_factor=1,
    h_spacing_factor = 1,
    **kwargs
):
    if 'show' in kwargs:
        plot_show = kwargs['show']
    else:
        plot_show = True
    kwargs['show'] = False
    if 'wspace' not in kwargs:
        kwargs['wspace'] = 0.25
    if 'hspace' not in kwargs:
        kwargs['hspace'] = 0.25

    n_rows = int(np.ceil(len(kwargs['keys']) / n_cols))
    fig_new, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4*n_cols*w_spacing_factor, 4*n_rows*h_spacing_factor)
    )
    fig_new.subplots_adjust(
        wspace=kwargs['wspace'],
        hspace=kwargs['hspace']
    )

    all_keys = kwargs['keys']
    del kwargs['keys']
    k = 0
    if n_rows == 1:
        for j in range(n_cols):
            if k<len(all_keys):
                sc.pl.violin(
                    adata,
                    keys=all_keys[k],
                    ax=axes[j],
                    **kwargs
                )
            else:
                axes[j].axis('off')
            k+=1
    elif n_cols == 1:
        for i in range(n_rows):
            if k<len(all_keys):
                sc.pl.violin(
                    adata,
                    keys=all_keys[k],
                    ax=axes[i],
                    **kwargs
                )
            else:
                axes[i].axis('off')
            k+=1
    else:
        for i in range(n_rows):
            for j in range(n_cols):
                if k<len(all_keys):
                    sc.pl.violin(
                        adata,
                        keys=all_keys[k],
                        ax=axes[i,j],
                        **kwargs
                    )
                else:
                    axes[i,j].axis('off')
                k+=1

    if plot_show: plt.show()

def __add_groupby_to_adata_gex_obs(adata, kwargs):
    if 'groupby' in kwargs:
        groupby = kwargs['groupby']
        if groupby not in adata.uns['gex_data'].obs.columns:
            if groupby in adata.obs.columns:
                adata.uns['gex_data'].obs[groupby] = adata.obs[groupby][
                    __get_indices_for_B_same_order_as_A(
                        adata.uns['gex_data'].obs_names.values,
                        adata.obs_names.values
                    )
                ]

def __add_groupby_to_adata_combo_obs(adata_combo, adata, kwargs):
    if 'groupby' in kwargs:
        groupby = kwargs['groupby']
        if groupby in adata.obs.columns:
            adata_combo.obs[groupby] = adata.obs[groupby]
        elif groupby in adata.uns['gex_data'].obs.columns:
            adata_combo.obs[groupby] = adata.uns['gex_data'].obs[groupby][
                __get_indices_for_B_same_order_as_A(
                    adata.obs_names.values,
                    adata.uns['gex_data'].obs_names.values
                )
            ]

def __get_adata_comb_w_kwargs(
    adata,
    **kwargs
):
    adata_combo, keys_combo = __get_adata_comb(
        adata,
        kwargs['keys'],
        basis = "X_pca"
    )
    kwargs['keys'] = keys_combo
    __add_groupby_to_adata_combo_obs(adata_combo, adata, kwargs)

    return adata_combo, kwargs

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------
def pca(adata,
        *,
        plot_pax=True,
        plot_gex=False,
        cmap_pax="RdBu_r",
        cmap_gex="viridis",
        cmap_obs="inferno",
        **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.pca.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot VIPER stored in adata on adata.obsm['X_pca'].
    plot_gex : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_pca'].
    cmap_pax : default: "RdBu_r"
        cmap to use for visualizing VIPER proteins.
    cmap_gex : default: "viridis"
        cmap to use for visualizing stored gExpr.
    cmap_obs : default: "inferno"
        cmap to use for visualizing stored numeric obs.
    **kwargs
        Arguments to provide to the sc.pl.pca function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    _plot(
        sc.pl.pca,
        "X_pca",
        adata,
        plot_pax,
        plot_gex,
        cmap_pax,
        cmap_gex,
        cmap_obs,
        kwargs
    )

def umap(adata,
         *,
         plot_pax=True,
         plot_gex=False,
         cmap_pax="RdBu_r",
         cmap_gex="viridis",
         cmap_obs="inferno",
         **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.umap.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot VIPER stored in adata on adata.obsm['X_umap'].
    plot_gex : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_umap'].
    cmap_pax : default: "RdBu_r"
        cmap to use for visualizing VIPER proteins.
    cmap_gex : default: "viridis"
        cmap to use for visualizing stored gExpr.
    cmap_obs : default: "inferno"
        cmap to use for visualizing stored numeric obs.
    **kwargs
        Arguments to provide to the sc.pl.pca function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    _plot(
        sc.pl.umap,
        "X_umap",
        adata,
        plot_pax,
        plot_gex,
        cmap_pax,
        cmap_gex,
        cmap_obs,
        kwargs
    )

def tsne(adata,
         *,
         plot_pax=True,
         plot_gex=False,
         cmap_pax="RdBu_r",
         cmap_gex="viridis",
         cmap_obs="inferno",
         **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.tsne.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot VIPER stored in adata on adata.obsm['X_tsne'].
    plot_gex : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_tsne'].
    cmap_pax : default: "RdBu_r"
        cmap to use for visualizing VIPER proteins.
    cmap_gex : default: "viridis"
        cmap to use for visualizing stored gExpr.
    cmap_obs : default: "inferno"
        cmap to use for visualizing stored numeric obs.
    **kwargs
        Arguments to provide to the sc.pl.tsne function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    _plot(
        sc.pl.tsne,
        "X_tsne",
        adata,
        plot_pax,
        plot_gex,
        cmap_pax,
        cmap_gex,
        cmap_obs,
        kwargs
    )

def diffmap(adata,
            *,
            plot_pax=True,
            plot_gex=False,
            cmap_pax="RdBu_r",
            cmap_gex="viridis",
            cmap_obs="inferno",
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.diffmap.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot VIPER stored in adata on adata.obsm['X_diffmap'].
    plot_gex : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_diffmap'].
    cmap_pax : default: "RdBu_r"
        cmap to use for visualizing VIPER proteins.
    cmap_gex : default: "viridis"
        cmap to use for visualizing stored gExpr.
    cmap_obs : default: "inferno"
        cmap to use for visualizing stored numeric obs.
    **kwargs
        Arguments to provide to the sc.pl.diffmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    _plot(sc.pl.diffmap,
          "X_diffmap",
          adata,
          plot_pax,
          plot_gex,
          cmap_pax,
          cmap_gex,
          cmap_obs,
          kwargs
    )

def draw_graph(adata,
               *,
               plot_pax=True,
               plot_gex=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.draw_graph.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot VIPER stored in adata on adata.obsm['X_draw_graph_fa'] or adata.obsm['X_draw_graph_fr'].
    plot_gex : default: False
        Plot gExpr stored in adata.uns['gex_data'] on adata.obsm['X_draw_graph_fa'] or adata.obsm['X_draw_graph_fr'].
    **kwargs
        Arguments to provide to the sc.pl.draw_graph function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax is True:
        pax_kwargs = __get_pax_params(kwargs, cmap_pax)
        pax_kwargs = __remove_missing_items(adata, pax_kwargs)
        sc.pl.draw_graph(adata, **pax_kwargs)

    if plot_gex is True:
        uns_data_slot = 'gex_data'
        adata_stored = adata.uns[uns_data_slot].copy()

        if "X_draw_graph_fa" in adata_stored.obsm.keys():
            adata_stored.obsm["X_draw_graph_fa"] = adata.obsm["X_draw_graph_fa"]
        if "X_draw_graph_fr" in adata_stored.obsm.keys():
            adata_stored.obsm["X_draw_graph_fr"] = adata.obsm["X_draw_graph_fr"]
        adata = adata_stored

        kwargs = __remove_missing_items(adata, kwargs)
        sc.pl.draw_graph(adata, **kwargs)

def spatial(adata,
            *,
            plot_pax=True,
            plot_gex=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.spatial.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata on adata.uns['spatial'].
    plot_gex : default: False
        Plot adata.uns['gex_data'] on adata.uns['spatial'].
    **kwargs
        Arguments to provide to the sc.pl.spatial function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __remove_missing_items(adata, kwargs)
        sc.pl.spatial(adata, **pax_kwargs)
    if(plot_gex is True):
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=None,
                                                       uns_slot='spatial')
        kwargs = __remove_missing_items(adata, kwargs)
        sc.pl.spatial(adata, **kwargs)


def embedding(adata,
              *,
              basis,
              plot_pax=True,
              plot_gex=False,
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
    plot_pax : default: True
        Plot adata on adata.obsm[basis].
    plot_gex : default: True
        Plot adata.uns['gex_data'] on adata.obsm[basis].
    **kwargs
        Arguments to provide to the sc.pl.embedding function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __remove_missing_items(adata, kwargs)
        sc.pl.embedding(adata, basis, **pax_kwargs)
    if plot_gex:
        kwargs = __remove_missing_items(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata,
            uns_data_slot='gex_data',
            obsm_slot=basis
        )
        sc.pl.embedding(adata, basis, **kwargs)

def embedding_density(adata,
                      *,
                      basis='umap',
                      plot_pax=True,
                      plot_gex=False,
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
    plot_pax : default: True
        Plot adata on adata.obsm[basis].
    plot_gex : default: False
        Plot adata.uns['gex_data'] on adata.obsm[basis].
    **kwargs
        Arguments to provide to the sc.pl.embedding_density function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        sc.pl.embedding_density(adata, basis, **kwargs)
    if plot_gex:
        adata = __get_stored_uns_data_and_prep_to_plot(adata,
                                                       uns_data_slot='gex_data',
                                                       obsm_slot=basis)
        sc.pl.embedding_density(adata, basis, **kwargs)

def scatter(adata,
            *,
            plot_pax=True,
            plot_gex=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.scatter.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.scatter function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __remove_missing_items(adata, kwargs)
        sc.pl.scatter(adata,**pax_kwargs)
    if plot_gex:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __remove_missing_items(adata, kwargs)
        sc.pl.scatter(adata,**kwargs)

def heatmap(adata,
            *,
            plot_pax=True,
            plot_gex=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.heatmap.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.heatmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __parse_var_names(adata, kwargs)
        sc.pl.heatmap(adata,**pax_kwargs)
    if plot_gex:
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.heatmap(adata,**kwargs)

def dotplot(adata,
            *,
            plot_pax=True,
            plot_gex=False,
            cmap_pax="Reds",
            cmap_gex="Greens",
            spacing_factor = 1,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.dotplot.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot VIPER stored in adata.
    plot_gex : default: False
        Plot gExpr stored in adata.uns['gex_data'].
    cmap_pax : default: "Reds"
        cmap to use for visualizing VIPER proteins.
    cmap_gex : default: "Greens"
        cmap to use for visualizing stored gExpr.
    spacing_factor : default: 1
        When plotting both pax and gex, adjust the size of the plot.
    **kwargs
        Arguments to provide to the sc.pl.dotplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax is True and plot_gex is True:
        _combo_dotplot(adata, cmap_pax, cmap_gex, spacing_factor, kwargs)
    elif plot_pax is True:
        pax_kwargs = kwargs.copy()
        if 'cmap' not in pax_kwargs:
            pax_kwargs['cmap'] = cmap_pax
        pax_kwargs = __parse_var_names(adata, pax_kwargs)
        sc.pl.dotplot(adata, **pax_kwargs)
    elif plot_gex is True:
        gex_kwargs = kwargs.copy()
        if 'cmap' not in gex_kwargs:
            gex_kwargs['cmap'] = cmap_gex
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        gex_kwargs = __parse_var_names(adata, gex_kwargs)
        sc.pl.dotplot(adata, **gex_kwargs)

def tracksplot(adata,
               *,
               plot_pax=True,
               plot_gex=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.tracksplot.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.tracksplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __parse_var_names(adata, kwargs)
        sc.pl.tracksplot(adata,**pax_kwargs)
    if plot_gex:
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.tracksplot(adata,**kwargs)

def violin(adata,
           *,
           plot_pax=True,
           plot_gex=False,
           n_cols = 4,
           w_spacing_factor=1,
           h_spacing_factor=1,
           **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.violin.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    n_cols : default: 4
        Number of columns in the violin plot.
    w_spacing_factor : default: 1
        Multiply the current width by this factor to stretch the plot.
    h_spacing_factor : default: 1
        Multiply the current height by this factor to stretch the plot.
    **kwargs
        Arguments to provide to the sc.pl.violin function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if isinstance(kwargs['keys'], str): kwargs['keys'] = [kwargs['keys']]

    if plot_pax and plot_gex:
        adata, kwargs = __get_adata_comb_w_kwargs(
            adata,
            **kwargs
        )
    elif plot_pax:
        kwargs = __parse_keys(adata, kwargs)
    elif plot_gex:
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_keys(adata, kwargs)

    _violin_grid(
        adata,
        n_cols,
        w_spacing_factor,
        h_spacing_factor,
        **kwargs
    )

def stacked_violin(adata,
                   *,
                   plot_pax=True,
                   plot_gex=False,
                   **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.stacked_violin.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.stacked_violin function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __parse_var_names(adata, pax_kwargs)
        sc.pl.stacked_violin(adata,**pax_kwargs)
    if plot_gex:
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.stacked_violin(adata,**kwargs)

def matrixplot(adata,
               *,
               plot_pax=True,
               plot_gex=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.matrixplot.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.matrixplot function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        pax_kwargs = __parse_var_names(adata, pax_kwargs)
        sc.pl.matrixplot(adata,**pax_kwargs)
    if plot_gex:
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        kwargs = __parse_var_names(adata, kwargs)
        sc.pl.matrixplot(adata,**kwargs)

def clustermap(adata,
               *,
               plot_pax=True,
               plot_gex=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.clustermap.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.clustermap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        sc.pl.clustermap(adata,**kwargs)
    if plot_gex:
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )
        sc.pl.clustermap(adata,**kwargs)

def ranking(adata,
            *,
            plot_pax=True,
            plot_gex=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.ranking.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.ranking function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        sc.pl.ranking(adata,**kwargs)
    if plot_gex:
        adata = __get_stored_uns_data_and_prep_to_plot(adata, uns_data_slot='gex_data')
        sc.pl.ranking(adata,**kwargs)

def dendrogram(adata,
               *,
               plot_pax=True,
               plot_gex=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.dendrogram.

    Parameters
    ----------
    adata
        Protein activity stored in an anndata object. Gene expression stored in
        adata.uns['gex_data'].
    plot_pax : default: True
        Plot adata.
    plot_gex : default: False
        Plot adata.uns['gex_data'].
    **kwargs
        Arguments to provide to the sc.pl.dendrogram function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if plot_pax:
        sc.pl.dendrogram(adata,**kwargs)
    if plot_gex:
        __add_groupby_to_adata_gex_obs(adata, kwargs)
        adata = __get_stored_uns_data_and_prep_to_plot(
            adata, uns_data_slot='gex_data'
        )

def ss_heatmap(
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
    _ss_heatmap(
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
    A wrapper around pyviper.pl.ss_heatmap to visualize the top master regulators (MRs) from each cluster using the VIPER anndata object.

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

def vis_net(
    net_Pruned,
    pax_data,
    gex_data,
    mr_list,
    cluster_labels = None,
    layout_alg='davidson_harel',
    seed = 0,
    figsize=(15, 15)
):
    """\
    Creates an igraph to visualize the relationship between regulators and targets of an Interactome object.

    Parameters
    ----------
    net_Pruned
        An object of class Interactome or a list of Interactome objects.
    pax_data
        Protein activity stored in an anndata object. Used to color the regulators.
    gex_data
        Gene expression stored in an anndata object. Used to color the targets.
    mr_list
        A list of proteins in net_Pruned to be visualized in the plot.
    cluster_labels : default: None
        Provide if you would like to create cluster-specific graphs.
    layout_alg : default: 'davidson_harel'
        Algorithm for the layout parameter from igraph.
    seed : default: 0
        Random seed for graph construction.
    figsize: default: (15, 15)
        figure dimensions for Matplotlib.
    """

    _vis_net(
        net_Pruned,
        pax_data,
        gex_data,
        mr_list,
        cluster_labels,
        layout_alg,
        seed,
        figsize
    )
