import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
import scanpy as sc
import pandas as pd
from random import sample
from random import seed
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from datetime import datetime

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ----------------------------- ** SIL DF FUNCS ** -----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def at_least_min(n, min):
    if n < min:
        return(min)
    else:
        return(n)
def update_sil_df(sil_df,
                  n_pcs,
                  a_res,
                  a_nn,
                  tot_clusters,
                  subsamp_iter,
                  silhouette_avgs_i,
                  seed,
                  curr_iter,
                  update_method = "iloc"):
    if(update_method == "iloc"):
        sil_df.iloc[curr_iter] = [curr_iter,
                                  n_pcs,
                                  a_res,
                                  a_nn,
                                  tot_clusters,
                                  subsamp_iter,
                                  silhouette_avgs_i,
                                  seed]
    elif(update_method == "concat"):
        new_iter_results = pd.DataFrame([[curr_iter,
                                          n_pcs,
                                          a_res,
                                          a_nn,
                                          tot_clusters,
                                          subsamp_iter,
                                          silhouette_avgs_i,
                                          seed]],
                                   columns=['iter',
                                            'n_pcs',
                                            'resolution',
                                            'knn',
                                            'n_clust',
                                            'subsamp_iter',
                                            'sil_avg',
                                            'seed'])
        sil_df = pd.concat([sil_df, new_iter_results],
                      ignore_index = True)
    return(sil_df)
def compute_silhouette_score(dist_mat, adata, pct_cells, SS_weights, SS_exp_base, i, obs_key_to_store_clusts = "clusters"):
    # In case Leiden only finds 1 cluster, set SS=0 so it won't be selected as optimal.
    if(len(adata.obs[obs_key_to_store_clusts].unique())==1):
        sil_avg = 0
    else:
        if SS_weights == 'exp':
            sil = silhouette_samples(X = dist_mat,
                                     labels = adata.obs[obs_key_to_store_clusts],
                                     metric="precomputed").tolist()
            if(pct_cells < 100):
                seed(i)
                sil = np.array(sample(sil, int(pct_cells/100*adata.shape[0])), dtype=float)
            else:
                sil = np.array(sil, dtype=float)

            neg_sil_indices = np.where(sil < 0)
            sil[neg_sil_indices] = -1 * (SS_exp_base ** np.abs(sil[neg_sil_indices]))
            sil_avg = np.mean(sil)
        else:
            sil_avg = silhouette_score(
                                            X = dist_mat,
                                            labels = adata.obs[obs_key_to_store_clusts],
                                            sample_size=int(pct_cells/100*adata.shape[0]), #Alessandro had commented out
                                            random_state = i,
                                            metric="precomputed"
        )
    return(sil_avg)
def add_clustering_results_to_sil_df_using_subsampling(dist_mat,
                                                       adata,
                                                       i,
                                                       pct_cells,
                                                       sil_df,
                                                       n_pcs,
                                                       a_res,
                                                       a_nn,
                                                       SS_weights,
                                                       SS_exp_base,
                                                       curr_iter,
                                                       update_method="iloc",
                                                       obs_key_to_store_clusts = "clusters"):
    silhouette_avgs_i = compute_silhouette_score(dist_mat, adata, pct_cells, SS_weights, SS_exp_base, i, obs_key_to_store_clusts)
    tot_clusters = len(adata.obs[obs_key_to_store_clusts].cat.categories)
    sil_df = update_sil_df(sil_df,
                           n_pcs,
                           a_res,
                           a_nn,
                           tot_clusters,
                           i,
                           silhouette_avgs_i,
                           i,
                           curr_iter,
                           update_method)
    return(sil_df)
def getEmptySilDF(nrow = 0):
    sil_df = pd.DataFrame(columns=['iter','n_pcs','resolution','knn','n_clust','subsamp_iter','sil_avg','seed'],
                         index=np.arange(nrow))
    return(sil_df)


# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ---------------------------- ** DISTANCE FUNCS ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def getDistanceMatrix(adata):
    from scipy.spatial.distance import pdist, squareform
    # make pairwise distance matrix
    pairwise_df = pd.DataFrame(
        squareform(pdist(adata, metric='euclidean')),
        columns = adata.index,
        index = adata.index
    )
    return(pairwise_df)

# SA_GS_subfunctions.R
def computeCorrDistance(data, verbose = True):
    # Equivalent to the following in R: d = sqrt(1 - stats::corr(X))
    if(verbose): print("Computing the correlation...")
    corr_matrix = np.corrcoef(data)
    if(verbose): print("Calculating sqrt_one_minus_corr_matrix...")
    sqrt_one_minus_corr_matrix = np.sqrt(1 - corr_matrix)
    if(verbose): print("Ensuring diagonal contains 0 values...")
    np.fill_diagonal(sqrt_one_minus_corr_matrix, 0)
    return(sqrt_one_minus_corr_matrix)
def add_row_column_names_to_dist_mat(dist_mat, adata):
    # Turn into a dataframe with row and column names
    df_dist = pd.DataFrame(
        dist_mat,
        columns = adata.obs_names,
        index = adata.obs_names
    )
    return(df_dist)
def getObjectDist(adata, object_dist, use_reduction, reduction_slot, verbose = True):
    if object_dist is None:
        if use_reduction == False:
            # use original features
            #
            d = computeCorrDistance(adata.X, verbose)
            numPCs = None
        elif use_reduction == True:
            # use principal components
            X = adata.obsm[reduction_slot]
            d = computeCorrDistance(X, verbose)
            numPCs = X.shape[1]
        else:
            raise ValueError("reduction must be logical.")
        d = add_row_column_names_to_dist_mat(d, adata)
    else:
        d = object_dist
        if use_reduction == True:
            numPCs = adata.obsm[reduction_slot].shape[1]
            cell_dims = 1 # cells are along rows
        del object_dist
    res = {"d": d, "numPCs": numPCs}
    return res

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ----------------------------- ** CLUSTER FUNCS ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def get_clust_func(clust_alg):
    if(str.lower(clust_alg) == "louvain"):
        clust_func = sc.tl.louvain
    else:
        clust_func = sc.tl.leiden
    return(clust_func)
def cluster_adata(adata,
                  my_random_seed,
                  res,
                  clust_alg,
                  obs_key_to_store_clusts = "clusters"):
                clust_alg_func = get_clust_func(clust_alg)
                clust_alg_func(adata,
                               random_state=my_random_seed,
                               resolution=res,
                               key_added=obs_key_to_store_clusts)
                return(adata)
