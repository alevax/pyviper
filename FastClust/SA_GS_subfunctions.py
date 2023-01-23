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
#     sil_df = sil_df.append(pd.DataFrame([[n_pcs,
#                                           a_res,
#                                           a_nn,
#                                           tot_clusters,
#                                           cur_iter,
#                                           silhouette_avgs_i,
#                                           seed]],
#                                        columns=['n_pcs',
#                                                 'resolution',
#                                                 'n_neighbors',
#                                                 'n_clust',
#                                                 'iter',
#                                                 'sil_avg',
#                                                 'seed']),ignore_index=True)
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
                                     labels = adata.obs[obs_key_to_store_clusts]).tolist()
            seed(i)
            sil = np.array(sample(sil, int(pct_cells/100*adata.shape[0])), dtype=float)
            neg_sil_indices = np.where(sil < 0)
            sil[neg_sil_indices] = -1 * (SS_exp_base ** np.abs(sil[neg_sil_indices]))
            sil_avg = np.mean(sil)
        else:
            sil_avg = silhouette_score(
                                            X = dist_mat,
                                            labels = adata.obs[obs_key_to_store_clusts],
                                            sample_size=int(pct_cells/100*adata.shape[0]), #Alessandro had commented out
                                            random_state = i
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
#     tmp = "Clusters: " + str(tot_clusters) + " | Avg Sil Score: " + str(round(silhouette_avgs_i,3))
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
#     distances = pdist(sample.values, metric='euclidean')
#     dist_matrix = squareform(distances)
#     return(dist_matrix)
    # make pairwise distance matrix
    pairwise_df = pd.DataFrame(
        squareform(pdist(adata, metric='euclidean')),
        columns = adata.index,
        index = adata.index
    )
    return(pairwise_df)

# SA_GS_subfunctions.R
# def cor_progbar(df): # RUNS TOO SLOW
#     # corr_matrix = np.zeros((len(df.columns), len(df.columns)))
#     # pbar = tqdm(desc="GridSearch",
#     #         total=len(df.columns) * len(df.columns),
#     #         position=0,
#     #         leave=True)
#     # for i in range(len(df.columns)):
#     #     for j in range(i, len(df.columns)):
#     #         corr_matrix[i, j] = df[df.columns[i]].corr(df[df.columns[j]])
#     #         pbar.update(1)
#     # pbar.close()
#     corr_matrix = np.zeros((len(df.columns), len(df.columns)))
#     with tqdm(total=len(df.columns)) as pbar:
#         for i in range(len(df.columns)):
#             corr_matrix[i] = df.corrwith(df[df.columns[i]])
#             pbar.update(1)
#     return(corr_matrix)
# def one_minus_x_progbar(x, pbar):
#     pbar.update()
#     return(1 - x)
# def sqrt_progbar(x, pbar):
#     pbar.update()
#     return(np.sqrt(x))
# def sqrt_one_minus_x_progbar(x, pbar):
#     pbar.update()
#     return(np.sqrt(1 - x))

def computeCorrDistance(data, verbose = True):
    # Equivalent to the following in R: d = sqrt(1 - stats::corr(X))
    # Convert numpy.array into a pandas dataframe
    # df = pd.DataFrame(X)
    # Compute the correlation (i.e. stats::corr(X) )
    # if(verbose): print("Computing the correlation...")
    # corr_matrix_df = df.corr() #cor_progbar(df)
    # Subtract every correlation value from 1 (i.e. 1 - corr_matrix )
    # if(verbose): print("Computing subtraction...")
    # with tqdm(total=len(corr_matrix_df.index) * len(corr_matrix_df.columns),
    #           unit='iteration',
    #           unit_scale=True,
    #           miniters=1,
    #           desc='Processing',
    #           leave=True) as pbar:
    #     df_sub = corr_matrix_df.applymap(lambda x: one_minus_x_progbar(x, pbar))
    # # Compute the square root of the previous result (i.e. sqrt(sub_matrix))
    # if(verbose): print("Computing square root...")
    # with tqdm(total=len(df_sub.index) * len(df_sub.columns),
    #           unit='iteration',
    #           unit_scale=True,
    #           miniters=1,
    #           desc='Processing',
    #           leave=True) as pbar:
    #     df_sqrt = df_sub.applymap(lambda x: sqrt_progbar(x, pbar))
    # if(verbose): print("Computing subtraction and square root...")
    # with tqdm(total=len(corr_matrix_df.index) * len(corr_matrix_df.columns),
    #           unit='iteration',
    #           unit_scale=True,
    #           miniters=1,
    #           desc='Processing',
    #           leave=True) as pbar:
    #     df_sqrt = corr_matrix_df.applymap(lambda x: sqrt_one_minus_x_progbar(x, pbar))
    # convert the dataframe to numpy array

    # In general, the performance of np.corrcoef will be faster when
        #  working with large datasets compared to df.corr() from pandas
    # if(verbose): print("Computing the correlation...")
    # with tqdm_numpy(desc="Calculating Correlation Matrix") as pbar:
    #     corr_matrix = np.corrcoef(X, rowvar=False)
    # if(verbose): print("Calculating sqrt_one_minus_corr_matrix...")
    # with tqdm_numpy(desc="Calculating sqrt_one_minus_corr_matrix") as pbar:
    #     sqrt_one_minus_corr_matrix = np.sqrt(1 - corr_matrix)

    # print("VERSION 1:")
    # corrcoef_start_time = time.time()
    # corr_matrix = np.corrcoef(data)
    # corrcoef_end_time = time.time()
    # corrcoef_elapsed_time = corrcoef_end_time - corrcoef_start_time
    # print(str(corrcoef_elapsed_time))
    # print("")

    # print("VERSION 2:")
    # corrcoef_start_time = time.time()
    # with tqdm(total=data.shape[1], desc="Calculating Correlation Matrix") as pbar:
    #     for i in range(data.shape[1]):
    #         corr_matrix = np.corrcoef(data, rowvar=False)
    #         pbar.update()
    # corrcoef_end_time = time.time()
    # corrcoef_elapsed_time = corrcoef_end_time - corrcoef_start_time
    # print(str(corrcoef_elapsed_time))
    if(verbose): print("Computing the correlation...")
    corr_matrix = np.corrcoef(data)
    if(verbose): print("Calculating sqrt_one_minus_corr_matrix...")
    sqrt_one_minus_corr_matrix = np.sqrt(1 - corr_matrix)
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
            d = computeCorrDistance(X, verbose) #computeCorrDistance(X.T, verbose)
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
