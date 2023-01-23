from SA_GS_subfunctions import *
from tqdm import tqdm
import seaborn as sns
from datetime import datetime

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# --------------------------- ** GRIDSEARCH FUNCS ** ---------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def get_gs_results(
    adata,
    res_vector=np.arange(0.1, 2, 0.2),
    NN_vector=np.arange(11, 101, 10),
    object_dist=None,
    use_reduction=True,
    reduction_slot="X_pca",
    SS_weights="unitary",
    SS_exp_base=2.718282,
    verbose=True,
    show_progress_bar=True,
    clust_alg="Leiden",
    n_subsamples=1,
    subsamples_pct_cells=100
):
    if verbose:
        print("Computing distance object...")
    dist_list = getObjectDist(adata, object_dist, use_reduction, reduction_slot)
    adata.obsm["dist_obj"] = dist_list["d"]
    n_pcs =  dist_list["numPCs"]
    n_iters = len(NN_vector)*len(res_vector)*n_subsamples
    sil_df = getEmptySilDF(n_iters)
    curr_iter = 0
    if verbose: print("Beginning GridSearch clustering...")
    if show_progress_bar: pbar = tqdm(desc = "GridSearch", total = n_iters, position=0, leave=True)
    for a_nn in NN_vector:
        # startTime_sc_pp_neighbors = datetime.now() # TIME TESTING
        # print(str(startTime_sc_pp_neighbors) + ": sc.pp.neighbors: a_nn=" + str(a_nn) + " - starting...") # TIME TESTING

        sc.pp.neighbors(adata, n_neighbors=a_nn, use_rep = "X_pca")

        # endTime_sc_pp_neighbors = datetime.now() # TIME TESTING
        # print(str(endTime_sc_pp_neighbors) + ": sc.pp.neighbors: a_nn=" + str(a_nn) + " - done.") # TIME TESTING

        # diffTime_sc_pp_neighbors = endTime_sc_pp_neighbors - startTime_sc_pp_neighbors # TIME TESTING
        # print(str(diffTime_sc_pp_neighbors.total_seconds())+ ": sc.pp.neighbors: a_nn=" + str(a_nn) + " - diffTime") # TIME TESTING

        for a_res in res_vector:
            # print("a_res=" + str(a_res))
            adata = cluster_adata(adata,
                                  0,#my_random_seed,
                                  a_res,
                                  clust_alg,
                                  obs_key_to_store_clusts = "GS_clusters")
            # run 100 times, change the seed
            silhouette_avgs = []
            for i in range(1,n_subsamples+1):
                sil_df = add_clustering_results_to_sil_df_using_subsampling(
                                   adata.obsm["dist_obj"],
                                   adata,
                                   i,
                                   subsamples_pct_cells,
                                   sil_df,
                                   n_pcs,
                                   a_res,
                                   a_nn,
                                   SS_weights,
                                   SS_exp_base,
                                   curr_iter,
                                   obs_key_to_store_clusts = "GS_clusters"
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

def get_opt_res_knn_from_gs_results(gs_results):
    search_df = gs_results["search_df"]
    max_sil_avg = np.max(search_df["sil_avg"])
    search_df_opt_row = search_df[search_df["sil_avg"] >= max_sil_avg].iloc[0]
    opt_res = search_df_opt_row["resolution"]
    opt_knn = int(search_df_opt_row["knn"])
    opt_values = {"opt_res": opt_res, "opt_knn": opt_knn}
    return(opt_values)

# -------------------------- ** MAIN RUN FUNCTION ** ---------------------------
def run_fastclust_GS_clustering(
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
    gs_results = get_gs_results(
        adata,
        res_vector,
        NN_vector,
        object_dist,
        use_reduction,
        reduction_slot,
        SS_weights,
        SS_exp_base,
        verbose,
        show_progress_bar,
        clust_alg,
        n_subsamples,
        subsamples_pct_cells
    )
    opt_values = get_opt_res_knn_from_gs_results(gs_results)
    opt_res = opt_values["opt_res"]
    opt_knn = opt_values["opt_knn"]
    gs_results["opt_result"] = [opt_values["opt_res"], opt_values["opt_knn"]]
    sc.pp.neighbors(adata, n_neighbors=opt_knn, use_rep = "X_pca")
    adata = cluster_adata(adata,
                          0,#my_random_seed,
                          opt_res,
                          clust_alg,
                          obs_key_to_store_clusts = "GS_clusters")
    adata.GS_results_dict = gs_results
    return(adata)

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ---------------------------- ** PLOTTING FUNCS ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def get_GS_search_plot(adata, plot_type = "sil_avg"):
    heatmap_table = pd.pivot_table(adata.GS_results_dict["search_df"],
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
    ax = sns.heatmap(heatmap_table,
                     cmap = color_map,
                     cbar_kws={'label': cbar_label},
                     linewidths=1,
                     linecolor='black')
    ax.invert_yaxis()
    # for _, spine in ax.spines.items():
        # spine.set_visible(True)
    plt.xlabel('Resolution')
    plt.ylabel('Nearest Neighbors')
    plt.close()
    return(fig)
