import igraph as ig
import matplotlib.pyplot as plt
import random
from matplotlib.colors import rgb2hex
import pandas as pd
import numpy as np

def _get_mean_pact_expr_per_clust(pax_data, gex_data, net_table, cluster_labels):
    # Get the average activity for regulators
    regulators = list(set(net_table['regulator']))
    data = pax_data.to_df().copy()
    data['cluster'] = cluster_labels
    pax_mean = data.groupby('cluster').mean()[regulators]

    # Get the average activity for targets
    targets = list(set(net_table['target']))
    data = gex_data.to_df().copy()
    data['cluster'] = cluster_labels
    gex_mean = data.groupby('cluster').mean()[targets]

    return pax_mean, gex_mean

def _map_to_colors(arr, color_palette = "RdBu_r"):
    arr[arr > 10] = 10
    arr[arr < -10] = -10

    # Normalize the array to the range [0, 1]
    #norm = (arr - arr.min()) / (arr.max() - arr.min())
    norm = (arr/20) + 0.5

    # Get the colormap
    cmap = plt.get_cmap(color_palette)

    # Map normalized values to the colormap
    colors = cmap(norm)

    colors = [rgb2hex(rgb) for rgb in colors]

    return colors

def _get_color_dicts(pax_mean, gex_mean):
    # Get regulator/target colors based on activity/expression
    pax_colors = _map_to_colors(pax_mean, color_palette = "RdBu_r")
    gex_colors = _map_to_colors(gex_mean, color_palette = "PuOr_r")
    # Label each color
    pax_colors_dict = dict(zip(pax_mean.index.values, pax_colors))
    gex_colors_dict = dict(zip(gex_mean.index.values, gex_colors))

    return pax_colors_dict, gex_colors_dict

def _assign_colors_to_nodes_using_dicts(g, pax_colors_dict, gex_colors_dict):
    # Assign colors to nodes
    v_colors = np.zeros(len(g.vs)).astype(str)

    for i in range(len(g.vs)):
        name = g.vs['name'][i]
        if name in list(pax_colors_dict.keys()):
            v_colors[i] = pax_colors_dict[name]
        else:
            v_colors[i] = gex_colors_dict[name]

    g.vs['color'] = v_colors

def _assign_colors_to_nodes(g, pax_mean, gex_mean):
    pax_colors_dict, gex_colors_dict = _get_color_dicts(
        pax_mean,
        gex_mean
    )
    _assign_colors_to_nodes_using_dicts(
        g,
        pax_colors_dict,
        gex_colors_dict
    )

def _get_graph_from_cluster_data(net_table, pax_mean_clust_i, gex_mean_clust_i, size_target=2):
    regulators = list(set(net_table['regulator']))
    targets = list(set(net_table['target']))

    edges = list(zip(net_table['regulator'], net_table['target']))
    g = ig.Graph(directed=True)
    g.add_vertices(list(set(regulators).union(targets)))  # Add all unique nodes
    g.add_edges(edges)  # Add edges

    # Assign node attributes
    # g.vs['color'] = ['red' if v['name'] in regulators else 'blue' for v in g.vs]

    # Define colors and sizes for nodes
    size_ratio = 15/3.5
    #size_target = 2
    g.vs['size'] = [size_target*size_ratio if v['name'] in regulators else size_target for v in g.vs]

    # Assign colors based on activity/expression of regulators/targets
    _assign_colors_to_nodes(g, pax_mean_clust_i, gex_mean_clust_i)

    # Label regulators but not targets
    vertex_label = [v["name"] if v['name'] in regulators else "" for v in g.vs]
    g.vs['label'] = vertex_label

    return g

def _plot_blobs(g, layout_alg, figsize=(15, 15)):
    # Generate layout and size correctly
    layout = g.layout(layout=layout_alg)
    layout.fit_into((0, 0, 100, 100))

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    ig.plot(
        g,
        target=ax,
        layout=layout,
        bbox=(800, 800),
        vertex_size=g.vs['size'],
        vertex_color=g.vs["color"],
        vertex_label=g.vs['label'],
        vertex_label_family = "sans",
        edge_width=0.1,
        edge_arrow_size=0,
        edge_color = "gray"
    )

def _vis_net(
    net_Pruned,
    pax_data,
    gex_data,
    mr_list,
    cluster_labels = None,
    layout_alg='davidson_harel',
    seed = 0,
    figsize=(15, 15),
    size_target=2
):
    random.seed(seed)
    # Filter net_table down to the regulators
    net_table = net_Pruned.net_table.copy()
    net_table = net_table[net_table['regulator'].isin(mr_list)]

    # Get mean activity/expression of each target/regulator per cluster
    if cluster_labels is None:
        cluster_labels = pd.Series(np.zeros(pax_data.shape[0]))
        cluster_labels.index = pax_data.obs_names

    pax_mean, gex_mean = _get_mean_pact_expr_per_clust(
        pax_data,
        gex_data,
        net_table,
        cluster_labels
    )

    # Targets may have a positive or negative relationship with regulators
    # Take the average of these relationships for visualization
    target_mean_mor = net_table.groupby('target')['mor'].mean().reset_index().set_index('target')
    gex_mean = gex_mean * target_mean_mor.loc[gex_mean.columns].values.flatten()

    # Iterate for each cluster
    unique_clusts = np.unique(cluster_labels)
    n_clusts = len(unique_clusts)
    for i in range(n_clusts):
        clust_i = unique_clusts[i]

        pax_mean_clust_i = pax_mean.loc[clust_i]
        gex_mean_clust_i = gex_mean.loc[clust_i]

        g = _get_graph_from_cluster_data(
            net_table, pax_mean_clust_i, gex_mean_clust_i,
            size_target=size_target
        )
        #_plot_blobs(g, layout_alg, seed, figsize)
        _plot_blobs(g, layout_alg, figsize)
