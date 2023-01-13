from SA_GS_subfunctions import *
from tqdm import tqdm
import seaborn as sns

def GridSearch(
    adata,
    res_vector = np.arange(0.1, 2, 0.2),
    NN_vector = np.arange(11, 101, 10),
    object_dist = None,
    use_reduction=True,
    reduction_slot="X_pca",
    SS_weights = "unitary",
    SS_exp_base = 2.718282,
    verbose = True,
    show_progress_bar = True,
    clust_alg = "Leiden",
    n_subsamples = 1,
    subsamples_pct_cells = 100
):
    if verbose: print("Computing distance object...")
    dist_list = getObjectDist(adata, object_dist, use_reduction, reduction_slot)
    adata.obsm["GS_obj"] = dist_list["d"]
    n_pcs =  dist_list["numPCs"]
    n_iters = len(NN_vector)*len(res_vector)*n_subsamples
    sil_df = getEmptySilDF(n_iters)
    curr_iter = 0
    if verbose: print("Beginning GridSearch clustering...")
    if show_progress_bar: pbar = tqdm(desc = "GridSearch", total = n_iters, position=0, leave=True)
    for a_nn in NN_vector:
        sc.pp.neighbors(adata, n_neighbors=a_nn, use_rep = "GS_obj")
        for a_res in res_vector:
            adata = cluster_adata(adata,
                             0,#my_random_seed,
                             a_res,
                             clust_alg)
            # run 100 times, change the seed
            silhouette_avgs = []
            for i in range(1,n_subsamples+1):
                sil_df = add_clustering_results_to_sil_df_using_subsampling(
                                   adata.obsm["GS_obj"],
                                   adata,
                                   i,
                                   subsamples_pct_cells,
                                   sil_df,
                                   n_pcs,
                                   a_res,
                                   a_nn,
                                   SS_weights,
                                   SS_exp_base,
                                   curr_iter
                )
                if show_progress_bar: pbar.update(1)
                curr_iter = curr_iter + 1
    if show_progress_bar: pbar.close()
    sil_df["resolution"] = np.around(sil_df["resolution"].astype(np.double),3)#prevent 0.3 being 0.300000000004

    run_params = {
        "res_vector": res_vector,
        "NN_vector": NN_vector,
        "use_reduction": use_reduction,
        "reduction_slot": reduction_slot,
        "SS_weights": SS_weights,
        "SS_exp_base": SS_exp_base,
        "clust_alg": clust_alg,
        "n_subsamples": n_subsamples,
        "subsamples_pct_cells": subsamples_pct_cells,
    }
    gs_results = {
        "search_df": sil_df,
        "run_params": run_params
    }
    return(gs_results)

def get_gs_search_plot(gs_results, plot_type = "sil_avg"):
    heatmap_table = pd.pivot_table(gs_results["search_df"],
                           values=plot_type,
                           index=['knn'],
                           columns=['resolution'], aggfunc=np.sum)
    if(plot_type == "sil_avg"):
        color_map = "YlOrRd"#"inferno"
    else: #plot_type = "n_clust"
        color_map = "YlGnBu"#"viridis"
    cbar_label = None
    if(plot_type == "sil_avg"):
        cbar_label = "Ave Sil Score"
    elif(plot_type == "n_clust"):
        cbar_label = "Clusters"

    fig = plt.figure()
    ax = sns.heatmap(heatmap_table, cmap = color_map, cbar_kws={'label': cbar_label})
    ax.invert_yaxis()
    plt.xlabel('Resolution')
    plt.ylabel('Nearest Neighbors')
    plt.close()
    return(fig)
