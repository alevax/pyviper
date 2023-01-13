from SA_GS_subfunctions import *
from dual_annealing_with_progress_bar import *
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import dual_annealing
import seaborn as sns

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ------------------------- ** DUAL ANNEALING FUNCS ** -------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def SA_Clustering(adata,
                  res_range = [0.01, 2],
                  NN_range = [3,30],
                  object_dist = None,
                  use_reduction=True,
                  reduction_slot="X_pca",
                  SS_weights = "unitary",
                  SS_exp_base = 2.718282,
                  verbose = True,
                  clust_alg = "Leiden",
                  n_subsamples = 1,
                  subsamples_pct_cells=100,
                  maxiter = 1000,
                  initial_temp = 5230,
                  restart_temp_ratio = 2e-5,
                  visit = 2.62,
                  accept = -5.0,
                  maxfun = 1e7,
                  seed = 0):
    # par_init = NULL,
    # control = NULL,
    # lq = 0,
    # rng_seeds = c(1234,0)):
    if verbose: print("Computing distance object...")
    dist_list = getObjectDist(adata, object_dist, use_reduction, reduction_slot)
    adata.obsm["GS_obj"] = dist_list["d"]
    n_pcs =  dist_list["numPCs"]
    # In order to use global inside a nested function, we have to declare a
        # variable with a global keyword inside a nested function
    global sil_df
    global curr_iter
    sil_df = getEmptySilDF(nrow = 0)
    curr_iter = 1
    bounds = [res_range, NN_range]

    def objective(v, adata, n_subsamples, subsamples_pct_cells, n_pcs, clust_alg, SS_weights, SS_exp_base):
        a_res, a_nn = v
        a_nn = int(np.floor(a_nn))
        global sil_df
        global curr_iter
        sc.pp.neighbors(adata, n_neighbors=a_nn, use_rep = "GS_obj")
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
                               curr_iter,
                               update_method = "concat"
            )
        curr_iter = curr_iter + 1
        sil_avg = np.mean(sil_df["sil_avg"].tail(n_subsamples).values) #Get the mean silouette score for these subsamples
        return(sil_avg*-1)

    if verbose: print("Beginning Simulated Annealing clustering...")
    # perform the dual annealing search
    # opt_result = dual_annealing(func = objective,
    opt_result = dual_annealing_with_progress_bar(func = objective,
                                bounds = bounds,
                                args = (adata,
                                        n_subsamples,
                                        subsamples_pct_cells,
                                        n_pcs,
                                        clust_alg,
                                        SS_weights,
                                        SS_exp_base),
                                maxiter = maxiter,
                                initial_temp = initial_temp,
                                restart_temp_ratio = restart_temp_ratio,
                                visit = visit,
                                accept = accept,
                                maxfun = maxfun,
                                seed = seed)
    run_params = {
        "res_range": res_range,
        "NN_range": NN_range,
        "use_reduction": use_reduction,
        "reduction_slot": reduction_slot,
        "SS_weights": SS_weights,
        "SS_exp_base": SS_exp_base,
        "clust_alg": clust_alg,
        "n_subsamples": n_subsamples,
        "subsamples_pct_cells": subsamples_pct_cells,
        "maxiter": maxiter,
        "initial_temp": initial_temp,
        "restart_temp_ratio": restart_temp_ratio,
        "visit": visit,
        "accept": accept,
        "maxfun": maxfun,
        "seed": seed
    }
    sa_results = {
        "search_df": sil_df,
        "opt_result": opt_result,
        "run_params": run_params
    }
    return(sa_results)

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# -------------------------- ** PLOTTING FUNCTIONS ** --------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def get_inferno_10():
    inferno_colors_10 = ["#000004FF",
                         "#1B0C42FF",
                         "#4B0C6BFF",
                         "#781C6DFF",
                         "#A52C60FF",
                         "#CF4446FF",
                         "#ED6925FF",
                         "#FB9A06FF",
                         "#F7D03CFF",
                         "#FCFFA4FF"]
    return(inferno_colors_10)
def get_viridis_10():
    virids_colors_10 = ["#440154FF",
                        "#482878FF",
                        "#3E4A89FF",
                        "#31688EFF",
                        "#26828EFF",
                        "#1F9E89FF",
                        "#35B779FF",
                        "#6DCD59FF",
                        "#B4DE2CFF",
                        "#FDE725FF"]
    return(virids_colors_10)
def get_YlOrRd_9():
    YlOrRd_9 = ["#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"]
    return(YlOrRd_9)

def get_BWR_3():
    BWR_3 = ["blue", "white", "red"]
    return(BWR_3)

def get_YlGnBu_9():
    YlGnBu_9 = ["#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"]
    return(YlGnBu_9)

def get_RdPu_9():
    RdPu_9 = ["#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"]
    return(RdPu_9)
def get_cmap_for_sa_search_plot(plot_type):
    if(plot_type == "sil_avg"):
        cmap = LinearSegmentedColormap.from_list("", get_YlOrRd_9())#["red","violet","blue"])
    elif(plot_type == "iter"):
        cmap = LinearSegmentedColormap.from_list("", get_YlGnBu_9())#["red","violet","blue"])
    elif(plot_type == "n_clust"):
        cmap = LinearSegmentedColormap.from_list("", get_RdPu_9())
    else:
        cmap = LinearSegmentedColormap.from_list("", get_inferno_10())
    return(cmap)
def get_cbar_label_for_sa_search_plot(plot_type):
    if(plot_type == "sil_avg"):
        cbar_label = "Ave Sil Score"
    elif(plot_type == "iter"):
        cbar_label = "Iteration"
    elif(plot_type == "n_clust"):
        cbar_label = "Clusters"
    else:
        cbar_label = None
    return(cbar_label)
def create_scatter_plot_for_sa_search_plot(ax, search_df, plot_type):
    cmap = get_cmap_for_sa_search_plot(plot_type)
    ax.scatter(search_df["resolution"], search_df["knn"], c=search_df[plot_type], cmap = cmap)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Nearest Neighbors')
def create_cbar_for_sa_search_plot(ax, plot_type):
    cbar_label = get_cbar_label_for_sa_search_plot(plot_type)
    PCM=ax.get_children()[0] #matplotlib.collections.PathCollection
    cbar = plt.colorbar(PCM, ax=ax)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(cbar_label, rotation=270)
def create_countour_layer_for_sa_search_plot(ax, search_df):
    countour_layer = sns.kdeplot(x=search_df["resolution"],
                                 y=search_df["knn"],
                                 ax = ax,
                                 clip = [ax.get_xlim(),ax.get_ylim()])
def get_sa_search_plot(sa_results, plot_type = "sil_avg", plot_density = True):
    # https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
    # if(plot_type == "sil_avg"):
    #     cmap = LinearSegmentedColormap.from_list("", get_YlOrRd_9())#["red","violet","blue"])
    #     cbar_label = "Ave Sil Score"
    # elif(plot_type == "iter"):
    #     cmap = LinearSegmentedColormap.from_list("", get_YlGnBu_9())#["red","violet","blue"])
    #     cbar_label = "Iteration"
    # elif(plot_type == "n_clust"):
    #     cbar_label = "Clusters"
    #     cmap = LinearSegmentedColormap.from_list("", get_RdPu_9())
    # else:
    #     cmap = LinearSegmentedColormap.from_list("", get_inferno_10())
    #     cbar_label = None
    # search_df = sa_results["search_df"]
    #
    # fig = plt.figure()
    # plt.scatter(search_df["resolution"], search_df["knn"], c=search_df[plot_type], cmap = cmap)
    # cbar = plt.colorbar()
    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.set_ylabel(cbar_label, rotation=270)
    # plt.xlabel('Resolution')
    # plt.ylabel('Nearest Neighbors')
    # plt.close()
    # return(fig)
    # https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
    search_df = sa_results["search_df"]
    fig, ax = plt.subplots()
    create_scatter_plot_for_sa_search_plot(ax, search_df, plot_type)
    create_cbar_for_sa_search_plot(ax, plot_type)
    if(plot_density == True):
        create_countour_layer_for_sa_search_plot(ax, search_df)
    return(fig)
