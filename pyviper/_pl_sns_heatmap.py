import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples
import seaborn as sns
import anndata
from ._tl import _find_top_mrs

# ----------------------------------
# -------- GET TOP MR FUNCS --------
# ----------------------------------

def _sort_top_mrs(adata, cluster_column, n_top_mrs, method="stouffer"):
    # Find top MRs for each cell type
    top_mrs_df = _find_top_mrs(
        adata,
        obs_column_name=cluster_column,
        N=n_top_mrs,
        method=method,
        return_as_df=True
    )

    # Ensure cluster_column is categorical
    if not isinstance(adata.obs[cluster_column].dtype, pd.CategoricalDtype):
        adata.obs[cluster_column] = adata.obs[cluster_column].astype('category')

    cell_types = adata.obs[cluster_column].cat.categories

    sorted_mrs = []

    for cell_type in cell_types:
        column_name = "mr_" + cell_type + "_scores"

        # Get MRs for this cell type
        cell_type_mrs = top_mrs_df[column_name].tolist()

        cell_indices = adata.obs[cluster_column] == cell_type

        # Calculate sum of values for each MR for this cell type
        mr_sums = []
        for mr in cell_type_mrs:
            if mr is not None:
                mr_sum = np.sum(adata[cell_indices, mr].X)
                mr_sums.append((mr, mr_sum))

        # Sort MRs by sum, from highest to lowest
        sorted_mrs_for_type = sorted(mr_sums, key=lambda x: x[1], reverse=True)
        sorted_mrs.extend([mr for mr, _ in sorted_mrs_for_type if mr is not None])

    return sorted_mrs

def _parse_top_mrs_list(
    top_mrs_list,
    adata,
    cluster_column,
    n_top_mrs,
    method,
    mr_col,
    h_clust_rows
):
    if top_mrs_list is None and mr_col is not None:
        if mr_col in adata.obs.columns:
            if adata.obs[mr_col].dtype == 'bool':
                top_mrs_list = adata.obs[adata.obs[mr_col] == True].index.tolist()
            else:
                top_mrs_list = adata.obs[adata.obs[mr_col] != 0].index.tolist()
        else:
            raise ValueError('top_mrs_list is None and mr_col not in adata.obs.columns.')
    else:
        top_mrs_list = _sort_top_mrs(adata, cluster_column, n_top_mrs, method)

    if h_clust_rows is True: top_mrs_list = np.unique(top_mrs_list)

    return top_mrs_list

# ----------------------------------
# -------- GET SS FUNCS --------
# ----------------------------------

def _get_ss(adata, cluster_column):
    silhouette_scores = silhouette_samples(adata.obsm['X_pca'], adata.obs[cluster_column])
    silhouette_scores = pd.Series(silhouette_scores, index=adata.obs.index)
    return silhouette_scores

def _get_sorted_indices(adata, cluster_column):
    if cluster_column is not None:
        # Calculate silhouette scores
        silhouette_scores = _get_ss(adata, cluster_column)

        # Sort samples by silhouette score within each cluster
        adata.obs['silhouette_score'] = silhouette_scores
        adata.obs['sort_order'] = \
            adata.obs.groupby(cluster_column, observed=True)['silhouette_score'].rank(ascending=False)
        sorted_indices = adata.obs.sort_values([cluster_column, 'sort_order']).index
    else:
        sorted_indices = adata.obs.index
    return sorted_indices

# ----------------------------------
# -------- COLOR PALETTES --------
# ----------------------------------

def get_vega_30():
    vega_10_dark = ['#1f77b4',
                   '#ff7f0e',
                   '#2ca02c',
                   '#d62728',
                   '#9467bd',
                   '#8c564b',
                   '#e377c2',
                   '#7f7f7f',
                   '#bcbd22',
                   '#17becf']
    vega_10_light = ['#aec7e8',
                    '#ffbb78',
                    '#98df8a',
                    '#ff9896',
                    '#c5b0d5',
                    '#c49c94',
                    '#f7b6d2',
                    '#c7c7c7',
                    '#dbdb8d',
                    '#9edae5']
    vega_10_medium = ["#1f497d", # (medium dark blue)
                      "#d2691e", # (medium dark orange)
                      "#228b22", # (medium dark green)
                      "#a52a2a", # (medium dark red/brown)
                      "#483d8b", # (medium dark purple)
                      "#7b7b7b", # (medium grey)
                      "#ffd700", # (medium yellow/gold)
                      "#008080", # (medium teal)
                      "#da70d6", # (medium pink/purple)
                      "#ffa07a"] # (medium salmon/orange))
    vega_30_complete = vega_10_dark + vega_10_light + vega_10_medium
    return vega_30_complete

def get_AI_60():
    AI_dark_set = [
    "#ff2d2d", #(Dark Red)
    "#ff8c00", #(Dark Orange)
    "#ffd700", #(Dark Yellow)
    "#808000", #(Olive)
    "#008000", #(Dark Green)
    "#008080", #(Teal)
    "#00cccc", #(Dark Cyan)
    "#000080", #(Navy)
    "#4b0082", #(Deep Purple)
    "#8b008b", #(Dark Magenta)
    "#9c6615", #(Brick)
    "#a0522d", #(Chocolate)
    "#654321", #(Dark Brown)
    "#708090", #(Slate Gray)
    "#555555", #(Dark Gray)
    "#4682b4", #(Steel Blue)
    "#8fbc8f", #(Dark Sea Green)
    "#a0522d", #(Sienna)
    "#9932cc", #(Dark Orchid)
    "#191970" #(Midnight Blue)
    ]
    AI_light_set = [
    "#ffc0cb", #(Pale Pink)
    "#ffdab9", #(Peach)
    "#fffacd", #(Lemon Chiffon)
    "#d7d7a0", #(Light Olive)
    "#98fb98", #(Pale Green)
    "#20b2aa", #(Light Sea Green)
    "#afeeee", #(Pale Turquoise)
    "#89cff0", #(Baby Blue)
    "#e6e6fa", #(Lavender)
    "#d8bfd8", #(Thistle)
    "#d2b48c", #(Tan)
    "#f4a460", #(Sandy Brown)
    "#ffe4c4", #(Bisque)
    "#d3d3d3", #(Light Gray)
    "#f5f5f5", #(White Smoke)
    "#87ceeb", #(Sky Blue)
    "#f5fffa", #(Mint Cream)
    "#eee8aa", #(Pale Goldenrod)
    "#fff0f5", #(Lavender Blush)
    "#f5f5dc" #(Beige)
    ]
    AI_medium_set = [
    "#ff0000", #(Red)
    "#ff7f50", #(Coral)
    "#ffff00", #(Yellow)
    "#00ff00", #(Lime)
    "#00ff7f", #(Spring Green)
    "#00ffff", #(Aqua)
    "#1e90ff", #(Dodger Blue)
    "#8a2be2", #(Blue Violet)
    "#da70d6", #(Orchid)
    "#ffa500", #(Orange)
    "#ffd700", #(Gold)
    "#f0e68c", #(Khaki)
    "#c0c0c0", #(Silver)
    "#dcdcdc", #(Gainsboro)
    "#696969", #(Dim Gray)
    "#6495ed", #(Cornflower Blue)
    "#66cdaa", #(Medium Aquamarine)
    "#9370db", #(Medium Purple)
    "#ff6347", #(Tomato)
    "#cd853f" #(Peru)
    ]
    AI_60_complete = AI_dark_set + AI_light_set + AI_medium_set
    return AI_60_complete

def _get_ith_discrete_pallete(i, n_colors):
    i = i % 7
    if n_colors > 60:
        return sns.color_palette("hls", n_colors)
    elif n_colors > 30:
        return sns.color_palette(get_AI_60(), n_colors)
    elif n_colors > 8 or i == 6:
        return sns.color_palette(get_vega_30(), n_colors)
    elif i == 0:
        return sns.color_palette("Set2", n_colors)
    elif i == 1:
        return sns.color_palette("Set3", n_colors)
    elif i == 2:
        return sns.color_palette("Dark2", n_colors)
    elif i == 3:
        return sns.color_palette("Accent", n_colors)
    elif i == 4:
        return sns.color_palette("Pastel1", n_colors)
    else:# i == 5:
        return sns.color_palette("Pastel2", n_colors)

def _get_ith_continuous_cmap(i):
    i = i % 6
    if i == 0:
        return plt.cm.Blues
    elif i == 1:
        return plt.cm.Oranges
    elif i == 2:
        return plt.cm.Greens
    elif i == 3:
        return plt.cm.Reds
    elif i == 4:
        return plt.cm.Purples
    elif i == 5:
        return plt.cm.Greys

# ----------------------------------
# -------- GET COLOR FUNCS ---------
# ----------------------------------

def _get_cluster_color_dict(adata, cluster_column):
    if cluster_column is None:
        return None
    elif not isinstance(adata.obs[cluster_column].dtype, pd.CategoricalDtype):
        adata.obs[cluster_column] = adata.obs[cluster_column].astype('category')

    cluster_categories = adata.obs[cluster_column].cat.categories

    n_colors = len(cluster_categories)
    if n_colors > 60:
        cluster_palette = sns.color_palette("hls", n_colors)
    elif n_colors > 30:
        cluster_palette = sns.color_palette(get_AI_60(), n_colors)
    else:
        cluster_palette = sns.color_palette(get_vega_30(), n_colors)

    cluster_color_dict = dict(zip(
        cluster_categories, cluster_palette
    ))
    return cluster_color_dict

def _get_annotation_color_dicts(adata, obs_metadata):
    if obs_metadata is None:
        return None

    annotation_color_dicts = {}
    for i, annotation in enumerate(obs_metadata):
        if annotation in adata.obs.columns:
            if isinstance(adata.obs[annotation].dtype, pd.CategoricalDtype) or adata.obs[annotation].dtype == 'object':
                # Categorical data
                adata.obs[annotation] = adata.obs[annotation].astype("category")
                unique_values = adata.obs[annotation].cat.categories
                annotation_colors = _get_ith_discrete_pallete(i, len(unique_values))
                annotation_color_dict = dict(zip(unique_values, annotation_colors))
                annotation_color_dicts[annotation] = annotation_color_dict
            else:
                # Continuous data
                annotation_color_dicts[annotation] = _get_ith_continuous_cmap(i)

    return annotation_color_dicts

def __add_cluster_to_col_colors(adata, col_colors, cluster_column, cluster_color_dict):
    if cluster_column is not None:
        col_colors[cluster_column] = adata.obs[cluster_column].astype(str).map(cluster_color_dict)

def __add_annotations_to_col_colors(adata, col_colors, obs_metadata, annotation_color_dicts):
    if obs_metadata is not None:
        for annotation in obs_metadata:
            if isinstance(adata.obs[annotation].dtype, pd.CategoricalDtype) or adata.obs[annotation].dtype == 'object':
                col_colors[annotation] = adata.obs[annotation].astype(str).map(annotation_color_dicts[annotation])
            else:
                cmap = annotation_color_dicts[annotation]
                # Normalize between 0 and 1 for plotting purposes
                min_val = adata.obs[annotation].min()
                max_val = adata.obs[annotation].max()
                norm = plt.Normalize(min_val, max_val)
                col_colors[annotation] = adata.obs[annotation].map(lambda x: cmap(norm(x)))

def __add_gexpr_to_col_colors(adata, col_colors, gexpr_metadata, plot_gex_norm):
    if gexpr_metadata is not None:
        if plot_gex_norm is True:
            gex_df = _get_gex_norm(adata)
        else:
            gex_df = adata.uns['gex_data'].to_df()
        for gene in gexpr_metadata:
            cmap = plt.cm.viridis
            # Normalize between 0 and 1 for plotting purposes
            min_val = gex_df[gene].min()
            max_val = gex_df[gene].max()
            norm = plt.Normalize(min_val, max_val)
            col_colors[gene + "_gexpr"] = gex_df[gene].map(lambda x: cmap(norm(x)))

def __add_viper_to_col_colors(adata, col_colors, viper_metadata):
    if viper_metadata is not None:
        pax_df = adata.to_df()
        for protein in viper_metadata:
            cmap = plt.cm.inferno
            # Normalize between 0 and 1 for plotting purposes
            min_val = pax_df[protein].min()
            max_val = pax_df[protein].max()
            norm = plt.Normalize(min_val, max_val)
            col_colors[protein] = pax_df[protein].map(lambda x: cmap(norm(x)))

def _get_col_colors(
    adata,
    cluster_column,
    cluster_color_dict,
    obs_metadata,
    annotation_color_dicts,
    gexpr_metadata,
    viper_metadata,
    plot_gex_norm
):
    col_colors = pd.DataFrame(index=adata.obs.index)

    # add cluster colors
    __add_cluster_to_col_colors(adata, col_colors, cluster_column, cluster_color_dict)

    # add annotation colors
    __add_annotations_to_col_colors(
        adata,
        col_colors,
        obs_metadata,
        annotation_color_dicts
    )

    # add gexpr colors
    __add_gexpr_to_col_colors(adata, col_colors, gexpr_metadata, plot_gex_norm)

    # add viper colors
    __add_viper_to_col_colors(adata, col_colors, viper_metadata)

    return col_colors

# ----------------------------------
# -------- MAIN HEATMAP FUNC ---------
# ----------------------------------

def _get_basic_clustermap(
    hm_data_sorted,
    var_names,
    sorted_indices,
    h_clust_rows,
    h_clust_cols,
    col_colors,
    show_gex_heatmap
):
    # Adjust fig size
    n_samps = hm_data_sorted.shape[1]
    n_vars = hm_data_sorted.shape[0]

    n_metadata = col_colors.shape[1]

    fig_width = n_samps*(10/8000) #10 #fig_width = 25
    fig_height = 0.5*(n_vars+n_metadata) + (0.3 if show_gex_heatmap else 0) #15 + (0.3 if show_gex_heatmap else 0)

    # colors_ratio is applied to each line of the cols_colors
    # In other words, if colors_ratio = 0.1, then each annotation in cols_colors
    # will be 10% of the total height.
    # prop_total_height_for_cols_colors = n_metadata/(n_metadata+n_vars)
    # Each col_colors should then have prop_total_height_for_cols_colors/n_metadata
    # so colors_ratio=1/(n_metadata+n_vars)

    cbar_y_pos = .6 if show_gex_heatmap else .4
    g = sns.clustermap(
        pd.DataFrame(
            hm_data_sorted,
            index=var_names,
            columns=sorted_indices
        ),
        row_cluster=h_clust_rows,
        col_cluster=h_clust_cols,
        col_colors=col_colors if n_metadata > 0 else None,
        colors_ratio=1/(n_metadata+n_vars),
        cmap="RdBu_r",
        center=0,
        figsize=(fig_width, fig_height),
        xticklabels=False,
        yticklabels=True,
        dendrogram_ratio=(0.15, 0.1),
        cbar_pos=(0.92, cbar_y_pos, 0.03, 0.35)
    )
    g.ax_heatmap.set_label('viper_clustermap')

    if h_clust_cols: g.ax_col_dendrogram.remove()

    # Adjust the position of the main heatmap
    _adjust_position_on_heatmap(g)

    return g

# ----------------------------------
# -------- PLOT HEATMAP FUNCS ---------
# ----------------------------------

def _adjust_position_on_heatmap(g):
    heatmap_pos = g.ax_heatmap.get_position()
    new_heatmap_pos = [heatmap_pos.x0, heatmap_pos.y0 - 0.02, heatmap_pos.width * 0.9, heatmap_pos.height]
    g.ax_heatmap.set_position(new_heatmap_pos)

def _get_gex_norm(adata):
    if adata.uns['gex_data'].raw is None:
        raise ValueError("adata.uns['gex_data'].raw is None. Need counts to plot normalized.")
    gex_raw = adata.uns['gex_data'].raw.to_adata()
    if np.all(gex_raw.X == np.floor(gex_raw.X)):
        sc.pp.normalize_total(gex_raw, target_sum=1e4)
        sc.pp.log1p(gex_raw)
    gex_data = gex_raw
    return gex_data

def _get_gex_heatmap_data(adata, var_names, sorted_indices, plot_gex_norm=False):
    if plot_gex_norm:
        gex_data = _get_gex_norm(adata)
    else:
        gex_data = adata.uns['gex_data']

    # Find the intersection of var_names with gex data, preserving order
    gex_genes = [gene for gene in var_names if gene in gex_data.var_names]

    # Prepare the gene expression data
    gex_data_subset = gex_data[:, gex_genes].X

    # Sort the columns (cells) to match the order in the main heatmap
    gex_data_sorted = gex_data_subset[gex_data.obs.index.get_indexer(sorted_indices), :]

    # Transpose the data for plotting (genes as rows, cells as columns)
    gex_data_sorted = gex_data_sorted.T

    return gex_data_sorted, gex_genes

def _add_gex_to_heatmap(
    g,
    adata,
    var_names,
    sorted_indices,
    plot_gex_norm,
    h_clust_rows
):
    # Reorder genes to be in same order as proteins
    if h_clust_rows:
        reordered_indices = g.dendrogram_row.reordered_ind
        var_names = [var_names[i] for i in reordered_indices]

    # Get the gene expression data using the final MR order
    gex_data_sorted, gex_genes = _get_gex_heatmap_data(
        adata,
        var_names,
        sorted_indices,
        plot_gex_norm
    )

    # Move other elements down
    for ax in g.fig.axes:
        bbox = ax.get_position()
        ax.set_position([bbox.x0, bbox.y0 - 0.3, bbox.width, bbox.height])

    heatmap_pos = g.ax_heatmap.get_position()

    # Add gene expression heatmap at the top
    ax_heatmap_gexpr = g.fig.add_axes([
        heatmap_pos.x0,
        heatmap_pos.y0 + heatmap_pos.height,
        heatmap_pos.width,
        heatmap_pos.height/len(var_names)*len(gex_genes)
        #heatmap_pos.height/len(var_names) gives us height per 1 item
        # and multiply by len(gex_genes) since may be different number
        # of genes than proteins
    ])#0.7
    sns.heatmap(gex_data_sorted,
                ax=ax_heatmap_gexpr,
                cmap="viridis",
                yticklabels=gex_genes,
                xticklabels=False,
                center=0,
                cbar_ax=g.fig.add_axes([.92, .73, .03, .3], label='gexpr_cbar'),
                )

    ax_heatmap_gexpr.set_label('gexpr_heatmap')

    ax_heatmap_gexpr.yaxis.set_label_position("right")
    ax_heatmap_gexpr.set_yticklabels(ax_heatmap_gexpr.get_yticklabels(), rotation=0)
    ax_heatmap_gexpr.yaxis.tick_right()


    gexpr_heatmap_pos = ax_heatmap_gexpr.get_position()

def _adjust_annotations(g, show_gex_heatmap):
    gexpr_heatmap_height = 0

    if show_gex_heatmap:
        gexpr_heatmap_ax = next(ax for ax in g.fig.axes if ax.get_label() == 'gexpr_heatmap')
        gexpr_heatmap_height = gexpr_heatmap_ax.get_position().height

    if g.ax_col_colors is not None:
        col_colors_pos = g.ax_col_colors.get_position()
        heatmap_pos = g.ax_heatmap.get_position()
        new_col_colors_pos = [
            heatmap_pos.x0,
            heatmap_pos.y0 + heatmap_pos.height + (gexpr_heatmap_height if show_gex_heatmap else 0),
            heatmap_pos.width,
            col_colors_pos.height
        ]
        g.ax_col_colors.set_position(new_col_colors_pos)

def _add_space_below_annotations(g):
    spacing = 0.005

    heatmap_pos = g.ax_heatmap.get_position()
    g.ax_heatmap.set_position([
        heatmap_pos.x0,
        heatmap_pos.y0 - spacing,
        heatmap_pos.width,
        heatmap_pos.height
    ])

    if g.ax_col_colors is not None:
        col_colors_pos = g.ax_col_colors.get_position()
        g.ax_col_colors.set_position([
            col_colors_pos.x0,
            col_colors_pos.y0 + spacing,
            col_colors_pos.width,
            col_colors_pos.height
        ])

def _adjust_row_dendrogram(g):
    row_dendrogram_pos = g.ax_row_dendrogram.get_position()
    heatmap_pos = g.ax_heatmap.get_position()
    g.ax_row_dendrogram.set_position([
        row_dendrogram_pos.x0,
        heatmap_pos.y0,
        row_dendrogram_pos.width,
        heatmap_pos.height
    ])

# ----------------------------------
# -------- ADD LEGENDS FUNCS ---------
# ----------------------------------

def _get_legend_vars_list(g):
    heatmap_pos = g.ax_heatmap.get_position()
    legend_width = 0.15
    legend_height = .225
    h_padding = 0.25 # h_padding = 0.05  # Horizontal padding
    legend_x = heatmap_pos.x0 + heatmap_pos.width + h_padding
    legend_y = heatmap_pos.y0 + heatmap_pos.height
    return [legend_x, legend_y, legend_width, legend_height]

# Add colorbars for continuous annotations
def _add_annotation_cbar(g, heatmap_pos, data, var_name, y_coord, cmap):
    legend_height = 0.02
    h_padding = 0.25
    legend_width = 0.15
    label_height = 0.075
    gap = 0.05

    cbar_x = heatmap_pos.x0 + heatmap_pos.width + h_padding
    cbar_ax = g.fig.add_axes([cbar_x, y_coord, legend_width, legend_height])
    norm = plt.Normalize(data[var_name].min(), data[var_name].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label(var_name + "_legend")
    return y_coord + legend_height + label_height + gap


def _add_annotation_legends(
    g,
    adata,
    cluster_column,
    cluster_color_dict,
    obs_metadata,
    annotation_color_dicts,
    gexpr_metadata,
    viper_metadata
):
    heatmap_pos = g.ax_heatmap.get_position()
    h_padding = 0.25
    legend_vars_list = _get_legend_vars_list(g)
    legend_x = legend_vars_list[0]
    legend_y = legend_vars_list[1]
    legend_width = legend_vars_list[2]
#     legend_height = legend_vars_list[3]

#     i = 0

    gap = 0.05
#     incr_factor = 0.025

#     ceiling = heatmap_pos.height#heatmap_pos.y0 + heatmap_pos.height

    y_coord = heatmap_pos.y0 + 0.1
    # Reverse order building up
    # VIPER annotations are last so add legends at bottom first
    viper_df = adata.to_df()
    if viper_metadata is not None:
        for protein in reversed(viper_metadata):
            y_coord = _add_annotation_cbar(
                g,
                heatmap_pos,
                data = viper_df,
                var_name = protein,
                y_coord = y_coord,
                cmap = plt.cm.inferno
            )

    # Stack gExpr annotation legends on top of VIPER annotations
    gexpr_df = adata.uns['gex_data'].to_df()
    if gexpr_metadata is not None:
        for gene in reversed(gexpr_metadata):
            y_coord = _add_annotation_cbar(
                g,
                heatmap_pos,
                data = gexpr_df,
                var_name = gene,
                y_coord = y_coord,
                cmap = plt.cm.viridis
            )

    y_coord = y_coord + gap

    # Stack obs_metadata annotations
    if obs_metadata is not None:
        # Calculate the number of discrete and continuous annotations
        discrete_annotations = [
            ann for ann in obs_metadata \
            if isinstance(adata.obs[ann].dtype, pd.CategoricalDtype) \
            or adata.obs[ann].dtype == 'object'
        ]
        continuous_annotations = [
            ann for ann in obs_metadata \
            if ann not in discrete_annotations
        ]

        # Stack continuous annotations
        for annotation in reversed(continuous_annotations):
            cbar_x = heatmap_pos.x0 + heatmap_pos.width + h_padding

            legend_height = 0.02
            cbar_ax = g.fig.add_axes([cbar_x, y_coord, legend_width, legend_height])
            norm = plt.Normalize(adata.obs[annotation].min(), adata.obs[annotation].max())
            sm = plt.cm.ScalarMappable(cmap=annotation_color_dicts[annotation], norm=norm)
            cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')

            label_obj = cbar.ax.text(0.5, 1.25, annotation, ha='center', va='bottom', transform=cbar.ax.transAxes)

            # Get the bounding box of the colorbar's Axes
            bbox_cbar = cbar_ax.get_window_extent().transformed(g.fig.transFigure.inverted())
            y_coord = y_coord + bbox_cbar.height*3 + gap


        for annotation in reversed(discrete_annotations):
            n_colors = len(annotation_color_dicts[annotation].values())
            legend_height = n_colors*0.01 + 0.01

            legend_x = heatmap_pos.x0 + heatmap_pos.width + h_padding
            legend_ax = g.fig.add_axes([legend_x, y_coord, legend_width, legend_height])
            handles = [plt.Rectangle((0,0),1,1, color=color)\
                        for color in annotation_color_dicts[annotation].values()]
            legend = legend_ax.legend(handles, annotation_color_dicts[annotation].keys(), title=annotation, loc='center')
            legend_ax.axis('off')
            legend_ax.set_label(annotation)

            bbox = legend.get_window_extent()
            bbox_normalized = bbox.transformed(g.fig.transFigure.inverted())

            legend_ax.set_position([legend_x, y_coord+bbox_normalized.height/2, legend_width, legend_height])

            y_coord = y_coord + bbox_normalized.height*1.5 + gap

    # Lastly, stack the cluster annotation on top
    if cluster_column is not None:
        n_colors = len(cluster_color_dict.values())
        legend_height = n_colors*0.01 + 0.01

        legend_ax = g.fig.add_axes([legend_x, y_coord, legend_width, legend_height])
        handles = [plt.Rectangle((0,0),1,1, color=color) for color in cluster_color_dict.values()]
        legend_ax.legend(handles, cluster_color_dict.keys(), title=cluster_column, loc='center')
        legend_ax.axis('off')
        legend_ax.set_label(cluster_column)

# –––––––––––––––––––––––––––––––––––––
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ########### MAIN FUNCTIONS ##########
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# –––––––––––––––––––––––––––––––––––––

def _general_heatmap(
    adata,
    var_names,
    cluster_column=None,
    show_gex_heatmap=False,
    plot_gex_norm=False,
    obs_metadata=None,
    gexpr_metadata=None,
    viper_metadata=None,
    h_clust_rows=True,
    h_clust_cols=False
):

    # Parse obs_metadata, gexpr_metadata, viper_metadata
    if isinstance(obs_metadata, str): obs_metadata = [obs_metadata]
    if isinstance(gexpr_metadata, str): gexpr_metadata = [gexpr_metadata]
    if isinstance(viper_metadata, str): viper_metadata = [viper_metadata]

    # Prepare the data for visualization
    hm_data = adata[:, var_names].X.T  # Extract and transpose data
    hm_data = np.array(hm_data, dtype=float)

    sorted_indices = _get_sorted_indices(adata, cluster_column)
    hm_data = hm_data[:, adata.obs.index.get_indexer(sorted_indices)]

    # Prepare color maps
    cluster_color_dict = _get_cluster_color_dict(adata, cluster_column)

    annotation_color_dicts = _get_annotation_color_dicts(adata, obs_metadata)
    col_colors = _get_col_colors(
        adata,
        cluster_column,
        cluster_color_dict,
        obs_metadata,
        annotation_color_dicts,
        gexpr_metadata,
        viper_metadata,
        plot_gex_norm
    )

    # Create the clustermap with seaborn
    g = _get_basic_clustermap(
        hm_data,
        var_names,
        sorted_indices,
        h_clust_rows = h_clust_rows,
        h_clust_cols = h_clust_cols,
        col_colors = col_colors,
        show_gex_heatmap = show_gex_heatmap
    )

    # Add gene expression heatmap if show_gex_heatmap is True
    if show_gex_heatmap:
        _add_gex_to_heatmap(
            g,
            adata,
            var_names,
            sorted_indices,
            plot_gex_norm,
            h_clust_rows
        )

    # Adjust the position of the metadata annotations
    _adjust_annotations(g, show_gex_heatmap)

    # Add spacing between annotation rows and data rows
    _add_space_below_annotations(g)

    # Adjust the position of the row dendrogram
    _adjust_row_dendrogram(g)

    # Add legends
    _add_annotation_legends(
        g,
        adata,
        cluster_column,
        cluster_color_dict,
        obs_metadata,
        annotation_color_dicts,
        gexpr_metadata,
        viper_metadata
    )

def _ss_heatmap(
    adata,
    var_names,
    cluster_column=None,
    show_gex_heatmap=False,
    plot_gex_norm=False,
    obs_metadata=None,
    gexpr_metadata=None,
    viper_metadata=None,
    h_clust_rows=True,
    h_clust_cols=False
):
    _general_heatmap(
        adata=adata,
        var_names=var_names,
        show_gex_heatmap=show_gex_heatmap,
        plot_gex_norm=plot_gex_norm,
        cluster_column=cluster_column,
        obs_metadata=obs_metadata,
        gexpr_metadata=gexpr_metadata,
        viper_metadata=viper_metadata,
        h_clust_rows=h_clust_rows,
        h_clust_cols=h_clust_cols
    )

def _mrs(
    adata,
    cluster_column,
    show_gex_heatmap=False,
    plot_gex_norm=False,
    obs_metadata=None,
    gexpr_metadata=None,
    viper_metadata=None,
    h_clust_rows=True,
    h_clust_cols=False,
    n_top_mrs=10,
    method='stouffer',
    top_mrs_list=None,
    mr_col=None
):
    # Get sorted top MRs
    if cluster_column is None:
        raise TypeError("cluster_column must be supplied to identify MRs.")

    top_mrs_list = _parse_top_mrs_list(
        top_mrs_list,
        adata,
        cluster_column,
        n_top_mrs,
        method,
        mr_col,
        h_clust_rows
    )

    if h_clust_rows:
        var_names = np.unique(top_mrs_list)
    else:
        var_names = top_mrs_list

    _general_heatmap(
        adata=adata,
        var_names=var_names,
        show_gex_heatmap=show_gex_heatmap,
        plot_gex_norm=plot_gex_norm,
        cluster_column=cluster_column,
        obs_metadata=obs_metadata,
        gexpr_metadata=gexpr_metadata,
        viper_metadata=viper_metadata,
        h_clust_rows=h_clust_rows,
        h_clust_cols=h_clust_cols
    )
