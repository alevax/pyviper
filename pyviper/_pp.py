import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from scipy.stats import rankdata, mannwhitneyu, norm
from scipy.special import ndtri, ndtr  #equivalent of norm.ppf and norm_cdf respectively but faster
import random
from tqdm import tqdm
from statsmodels.stats import multitest
from ._filtering_funcs import _get_anndata_filtered_by_feature_group
from tqdm.auto import tqdm
import os

def norm_ppf(x):
    return ndtri(x)

# def norm_cdf(x):
#     return ndtr(x)

def norm_sf(x):
    return ndtr(-x)

def _median(a, axis=None, out=None, overwrite_input=False, keepdims=False):
    return np.median(a, axis, out, overwrite_input, keepdims)

def _mad_from_R(x, center=None, constant=1.4826, low=False, high=False):
    if center is None:
        center=np.median(x)
    x = x[~np.isnan(x)] if np.isnan(x).any() else x
    n = len(x)
    if (low or high) and n % 2 == 0:
        if low and high:
            raise ValueError("'low' and 'high' cannot be both True")
        n2 = n // 2 + int(high)
        return constant * np.sort(np.abs(x - center))[n2]
    return constant * np.median(np.abs(x - center))

# Function assumes features as rows and observations as columns
# Numerator Functions:
    # Median - numpy.median
    # Mean - numpy.mean
# Denominator Functions:
    # Median absolute deviation - mad_from_R
    # Standard deviation - statistics.stdev
def _rank_norm_df(x, NUM_FUN=np.median, DEM_FUN = _mad_from_R, verbose = False):
    rank = rankdata(x, axis=0)
    median = NUM_FUN(rank, axis=1, keepdims=True)#np.median(rank, axis=1, keepdims=True)
    mad = np.apply_along_axis(DEM_FUN, 1, rank)

    x = ((rank - median)/mad[:, np.newaxis])

    if verbose: print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
    if verbose: print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
    if verbose: print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))
    if verbose: print("- Features to Remove:")

    # Take care of infinite values
    max_finite = np.nanmax(x[np.isfinite(x)])
    min_finite = np.nanmin(x[np.isfinite(x)])
    x[np.isposinf(x)] = max_finite
    x[np.isneginf(x)] = min_finite

    x = np.where(np.isnan(x), np.nanmin(x), x)
    x = np.clip(x, a_min=np.nanmin(x), a_max=np.nanmax(x))
    if verbose: print("- Removing NULL/NA features ...")
    x = x[~np.isnan(x).any(axis=1)]

    if verbose: print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
    if verbose: print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
    if verbose: print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))

    return x

def _rank_norm(
    adata,
    NUM_FUN = np.median,
    DEM_FUN = _mad_from_R,
    layer = None,
    key_added = None,
    copy = False
):
    if copy is True: adata = adata.copy()

    if isinstance(adata, pd.DataFrame) or isinstance(adata, np.ndarray):
        adata[:] = _rank_norm_df(adata, NUM_FUN, DEM_FUN)
    elif(isinstance(adata, anndata.AnnData) or isinstance(adata, anndata._core.anndata.AnnData)):
        if layer is None:
            input_array = adata.X
        else:
            input_array = adata.layers[layer]

        transformed_array = _rank_norm_df(input_array, NUM_FUN, DEM_FUN)

        if key_added is not None:
            adata.layers[key_added] = transformed_array
        elif layer is not None:
            adata.layers[layer] = transformed_array
        else:
            adata.X = transformed_array
    else:
        raise Exception("In RankNorm(x), x must be anndata.AnnData, numpy.ndarray or pandas.DataFrame.")

    if copy is True:
        return adata
    else:
        return None

def __sigT(x, slope = 20, inflection = 0.5):
    return (1 - 1/(1 + np.exp(slope * (x - inflection))))









# Interleave the columns
def _interleave_dfs(df_a, df_b):
    # Interleave the rows
    interleaved_rows = [i for pair in zip(df_a.index, df_b.index) for i in pair]
    merge_df = pd.concat([df_a, df_b], axis=0, ignore_index=False)#ignore_index=True)
    merge_df = merge_df.loc[interleaved_rows]
    return merge_df

def _corr_btw_2D_and_1D_arrays(arr_2d, arr_1d, axis):
    # Subtract mean along the specified axis
    arr_2d_mean = np.mean(arr_2d, axis=axis, keepdims=True)
    arr_1d_mean = np.mean(arr_1d)
    arr_2d_centered = arr_2d - arr_2d_mean
    arr_1d_centered = arr_1d - arr_1d_mean

    # Compute standard deviations
    arr_2d_std = np.std(arr_2d, axis=axis, keepdims=True)
    arr_1d_std = np.std(arr_1d)

    # Compute correlation coefficient
    if axis == 0:
        corr = np.sum(arr_2d_centered * arr_1d_centered[:, np.newaxis], axis=axis) / (arr_2d.shape[0] - 1)
        corr /= (arr_2d_std.squeeze() * arr_1d_std)
    elif axis == 1:
        corr = np.sum(arr_2d_centered * arr_1d_centered, axis=axis) / (arr_2d.shape[1] - 1)
        corr /= (arr_2d_std * arr_1d_std).squeeze()
    else:
        raise ValueError("axis must be 0 or 1")

    return corr

def _get_closeness_to_centroid(pca_array, cluster_mask):
    clust_centroid = np.mean(pca_array[cluster_mask], axis = 0)
    corr_values = _corr_btw_2D_and_1D_arrays(arr_2d = pca_array, arr_1d = clust_centroid, axis = 1)
    sqrt_one_minus_corr_values = np.sqrt(1 - corr_values)
    distances = sqrt_one_minus_corr_values
    closeness = 1 - distances
    return closeness

def _get_var_corr_with_clust(dat_array_ranked, pca_array, cluster_mask):
    closeness_to_clust_centroid = _get_closeness_to_centroid(pca_array, cluster_mask)
    corr_with_clust = _corr_btw_2D_and_1D_arrays(
        arr_2d = dat_array_ranked,
        arr_1d = closeness_to_clust_centroid,
        axis = 0
    )
    return corr_with_clust

def _corr_btw_2D_and_1D_arrays_axis0_precomp2D(arr_1d, arr_2d, arr_2d_centered, arr_2d_std_sqz):
    # Subtract mean along the specified axis
    arr_1d_mean = np.mean(arr_1d)
    arr_1d_centered = arr_1d - arr_1d_mean
    # Compute standard deviations
    arr_1d_std = np.std(arr_1d)

    # Compute correlation coefficient
    corr = np.sum(arr_2d_centered * arr_1d_centered[:, np.newaxis], axis=0) / (arr_2d.shape[0] - 1)
    corr /= (arr_2d_std_sqz * arr_1d_std)
    return corr

def _spearman_clusters_df(dat_df, pca_array, cluster_vector, compute_pvals = True, null_iters = 1000, verbose = True):
    # Convert the DataFrame to a NumPy array
    dat_array = dat_df.to_numpy()
    dat_array_ranked = rankdata(dat_array, axis = 0)

    # Find unique clusters and initialize arrays to store Spearman scores
    unique_clusters, cluster_indices = np.unique(cluster_vector, return_inverse=True)
    n_clusters = len(unique_clusters)
    n_genes = dat_df.shape[1]
    spearman_scores_arr = np.zeros((n_clusters, n_genes))
    if compute_pvals: spearman_pvals_arr = np.zeros((n_clusters, n_genes))

    # Calculate Spearman scores for each cluster and gene
    for i in tqdm(range(n_clusters)) if verbose else range(n_clusters):
        cluster_mask = (cluster_indices == i)

        spearman_scores = _get_var_corr_with_clust(dat_array_ranked, pca_array, cluster_mask)
        spearman_scores_arr[i, :] = spearman_scores

        if compute_pvals:
            # Compute null model
            spearman_scores_null_dist = np.zeros((null_iters, n_genes))
            closeness_to_clust_centroid = _get_closeness_to_centroid(pca_array, cluster_mask)

            # Precompute values to speed up the computation of the Spearman
                # Subtract mean along the specified axis
            dat_array_ranked_mean = np.mean(dat_array_ranked, axis=0, keepdims=True)
            dat_array_ranked_centered = dat_array_ranked - dat_array_ranked_mean
                # Compute standard deviations
            dat_array_ranked_std = np.std(dat_array_ranked, axis=0, keepdims=True)

            for j in range(null_iters):
                random.seed(j)
                closeness_to_clust_centroid_rand = random.sample(
                    list(closeness_to_clust_centroid),
                    len(closeness_to_clust_centroid)
                )

                spearman_scores_null_dist[j, :] = _corr_btw_2D_and_1D_arrays_axis0_precomp2D(
                    arr_1d = closeness_to_clust_centroid_rand,
                    arr_2d = dat_array_ranked,
                    arr_2d_centered = dat_array_ranked_centered,
                    arr_2d_std_sqz = dat_array_ranked_std.squeeze()
                )

            mean_ests = np.mean(spearman_scores_null_dist, axis=0)
            sd_ests = np.std(spearman_scores_null_dist, axis=0)
            pvals = 2 * norm.cdf(-abs(spearman_scores), loc=mean_ests, scale=sd_ests)
            pvals[pvals > 1] = 1
            spearman_pvals_arr[i, :] = pvals

    # Create a DataFrame from the computed Stouffer scores
    spearman_scores_df = pd.DataFrame(spearman_scores_arr, index=unique_clusters + "_scores", columns=dat_df.columns)
    if compute_pvals:
        spearman_pvals_df = pd.DataFrame(spearman_pvals_arr, index=unique_clusters + "_pvals", columns=dat_df.columns)
        result_df = _interleave_dfs(spearman_scores_df, spearman_pvals_df)
    else:
        result_df = spearman_scores_df

    return result_df

def _stouffer_clusters_df(dat_df, cluster_vector, compute_pvals = True, null_iters = 1000, verbose = True):
    return _gen_sig_clusters_df(dat_df, cluster_vector, compute_pvals, null_iters, verbose, method = "stouffer")

def _mean_diffs_clusters_df(dat_df, cluster_vector, compute_pvals = True, null_iters = 1000, verbose = True):
    return _gen_sig_clusters_df(dat_df, cluster_vector, compute_pvals, null_iters, verbose, method = "mean_diffs")

def _gen_sig_clusters_df(dat_df, cluster_vector, compute_pvals = True, null_iters = 1000, verbose = True, method = "stouffer"):

    # Convert the DataFrame to a NumPy array
    dat_array = dat_df.to_numpy()

    # Find unique clusters and initialize arrays to store Stouffer scores
    unique_clusters, cluster_indices = np.unique(cluster_vector, return_inverse=True)
    n_clusters = len(unique_clusters)
    n_genes = dat_df.shape[1]
    scores_arr = np.zeros((n_clusters, n_genes))
    if compute_pvals: pvals_arr = np.zeros((n_clusters, n_genes))

    # Calculate the denominator for Stouffer scores for each cluster
    cluster_sizes = np.bincount(cluster_indices)
    sqrt_cluster_sizes = np.sqrt(cluster_sizes)

    # Calculate Stouffer scores for each cluster and gene
    for i in tqdm(range(n_clusters)) if verbose else range(n_clusters):
        cluster_mask = (cluster_indices == i)
        cluster_data = dat_array[cluster_mask]
        other_data = dat_array[~cluster_mask]
        if method == "stouffer":
            scores = np.sum(cluster_data, axis=0) / sqrt_cluster_sizes[i]
        elif method == "mean_diffs":
            scores = np.mean(cluster_data, axis=0) - np.mean(other_data, axis=0)
        elif method == "median_diffs":
            scores = np.median(cluster_data, axis=0) - np.median(other_data, axis=0)

        scores_arr[i, :] = scores

        if compute_pvals:
            scores_null_dist = np.zeros((null_iters, n_genes))
            for j in range(null_iters):
                random.seed(j)
                cluster_mask_rand = np.array(random.sample(list(cluster_mask), len(cluster_mask)))
                cluster_data_rand = dat_array[cluster_mask_rand]
                other_data_rand = dat_array[~cluster_mask_rand]
                if method == "stouffer":
                    scores_null_dist[j, :] = np.sum(cluster_data_rand, axis=0) / sqrt_cluster_sizes[i]
                elif method == "mean_diffs":
                    scores_null_dist[j, :] = np.mean(cluster_data_rand, axis=0) - np.mean(other_data_rand, axis=0)
                elif method == "median_diffs":
                    scores_null_dist[j, :] = np.median(cluster_data_rand, axis=0) - np.median(other_data_rand, axis=0)

            mean_ests = np.mean(scores_null_dist, axis=0)
            sd_ests = np.std(scores_null_dist, axis=0)
            pvals = 2 * norm.cdf(-abs(scores), loc=mean_ests, scale=sd_ests)
            pvals[pvals > 1] = 1
            pvals_arr[i, :] = pvals

    # Create a DataFrame from the computed Stouffer scores
    scores_df = pd.DataFrame(scores_arr,
                             index=np.array([x + "_scores" for x in unique_clusters]),
                             columns=dat_df.columns)

    if compute_pvals:
        pvals_df = pd.DataFrame(pvals_arr,
                                index=np.array([x + "_pvals" for x in unique_clusters]),
                                columns=dat_df.columns)
        result_df = _interleave_dfs(scores_df, pvals_df)
    else:
        result_df = scores_df

    return result_df

def _mwu_clusters_df(dat_df, cluster_vector, compute_pvals = True, verbose = True):
    # Ensure cluster_vector has the same number of samples as rows in dat_df
    if len(cluster_vector) != dat_df.shape[0]:
        raise ValueError("Cluster vector length does not match the number of rows in the DataFrame.")
    unique_clusters, cluster_indices = np.unique(cluster_vector, return_inverse=True)
    n_clusters = len(unique_clusters)

    # Need at least 2 clusters to perform MWU
    if n_clusters < 2:
        raise ValueError("Only ", str(n_clusters), " unique clusters. At least 2 needed to compute MWU.")

    # Convert the DataFrame to a NumPy array
    dat_array = dat_df.to_numpy()

    # Find unique clusters and initialize arrays to store Stouffer scores
    unique_clusters, cluster_indices = np.unique(cluster_vector, return_inverse=True)

    n_genes = dat_df.shape[1]
    mwu_scores = np.zeros((n_clusters, n_genes))
    mwu_pvals = np.zeros((n_clusters, n_genes))

    # Calculate Stouffer scores for each cluster and gene
    for i in tqdm(range(n_clusters)) if verbose else range(n_clusters):
        clust_i_name = unique_clusters[i]
        cluster_mask = (cluster_vector == clust_i_name)
        clust_i_samples = dat_array[cluster_mask,]
        clust_others_samples = dat_array[~cluster_mask,]
        mwu_output = mannwhitneyu(clust_i_samples, clust_others_samples)
        mwu_scores[i, :] = mwu_output.statistic
        mwu_pvals[i, :] = mwu_output.pvalue

    # Create a DataFrame from the computed Stouffer scores
    scores_df = pd.DataFrame(mwu_scores, index=unique_clusters + "_scores", columns=dat_df.columns)
    pvals_df = pd.DataFrame(mwu_pvals, index=unique_clusters + "_pvals", columns=dat_df.columns)

    if compute_pvals:
        result_df = _interleave_dfs(scores_df, pvals_df)
    else:
        result_df = scores_df

    return result_df

def _sig_clusters_adata(adata,
                        obs_column_name = None,
                        layer = None,
                        filter_by_feature_groups = None,
                        sig_method = "stouffer",
                        compute_pvals = True,
                        null_iters = 1000,
                        verbose = True,
                        pca_slot = "X_pca"):
    if isinstance(adata, pd.DataFrame): adata = anndata.AnnData(adata)

    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    if layer is None:
        dat_df = adata_filt.to_df()
    else:
        dat_df = pd.DataFrame(adata_filt.layers[layer],
                              index=adata_filt.obs_names,
                              columns=adata_filt.var_names)
    if obs_column_name is None:
        cluster_vector = np.zeros(adata_filt.shape[0]).astype(int).astype(str)
    elif isinstance(obs_column_name, str):
        cluster_vector = adata.obs[obs_column_name]
    else:
        cluster_vector = obs_column_name

    if sig_method == "stouffer":
        result_df = _stouffer_clusters_df(dat_df, cluster_vector, compute_pvals, null_iters, verbose)
    elif sig_method == "mwu":
        result_df = _mwu_clusters_df(dat_df, cluster_vector, compute_pvals, verbose)
    elif sig_method == "spearman":
        result_df = _spearman_clusters_df(dat_df, adata.obsm[pca_slot], cluster_vector, compute_pvals, null_iters, verbose)
    elif sig_method == "mean_diffs":
        result_df = _mean_diffs_clusters_df(dat_df, cluster_vector, compute_pvals, null_iters, verbose)
    else:
        raise ValueError("sig_method must be stouffer, mwu, or spearman")

    return result_df

def _sig_clusters(adata,
                  obs_column_name = None,
                  layer = None,
                  filter_by_feature_groups = None,
                  key_added = 'stouffer',
                  return_as_df = False,
                  copy = False,
                  sig_method = 'stouffer',
                  compute_pvals = True,
                  null_iters = 1000,
                  verbose = True,
                  pca_slot = "X_pca"):
    if not isinstance(adata, anndata.AnnData):
        if isinstance(adata, pd.DataFrame):
            return_as_df = True
        else:
            raise ValueError("adata must be anndata.AnnData or pd.DataFrame.")
    elif copy is True and return_as_df is True:
        raise ValueError("copy and return_as_df cannot both be True when adata is anndata.AnnData.")

    if copy is True: adata = adata.copy()

    result_df = _sig_clusters_adata(adata,
                                    obs_column_name,
                                    layer,
                                    filter_by_feature_groups,
                                    sig_method,
                                    compute_pvals,
                                    null_iters,
                                    verbose,
                                    pca_slot)
    if return_as_df is True:
        return result_df
    else:
        result_df = result_df.T
        result_df.columns = key_added + "_" + result_df.columns
        adata.var = pd.concat([adata.var, result_df], axis=1, join='outer')
    if copy is True:
        return adata

def _stouffer(adata,
              obs_column_name = None,
              layer = None,
              filter_by_feature_groups = None,
              key_added = 'stouffer',
              compute_pvals = True,
              null_iters = 1000,
              verbose = True,
              return_as_df = False,
              copy = False):
    return _sig_clusters(
        adata,
        obs_column_name,
        layer,
        filter_by_feature_groups,
        key_added,
        return_as_df,
        copy,
        sig_method = 'stouffer',
        compute_pvals = compute_pvals,
        null_iters = null_iters,
        verbose = verbose,
        pca_slot = None
    )

def _mwu(adata,
         obs_column_name = None,
         layer = None,
         filter_by_feature_groups = None,
         key_added = 'mwu',
         compute_pvals = True,
         verbose = True,
         return_as_df = False,
         copy = False):
    return _sig_clusters(
        adata,
        obs_column_name,
        layer,
        filter_by_feature_groups,
        key_added,
        return_as_df,
        copy,
        sig_method = 'mwu',
        compute_pvals = compute_pvals,
        null_iters = None,
        verbose = verbose,
        pca_slot = None
    )

def _spearman(adata,
              pca_slot = "X_pca",
              obs_column_name = None,
              layer = None,
              filter_by_feature_groups = None,
              key_added = 'spearman',
              compute_pvals = True,
              null_iters = 1000,
              verbose = True,
              return_as_df = False,
              copy = False):
    return _sig_clusters(
        adata,
        obs_column_name,
        layer,
        filter_by_feature_groups,
        key_added,
        return_as_df,
        copy,
        sig_method = 'spearman',
        compute_pvals = compute_pvals,
        null_iters = null_iters,
        verbose = verbose,
        pca_slot = pca_slot
    )

def _mean_diffs(adata,
                obs_column_name = None,
                layer = None,
                filter_by_feature_groups = None,
                key_added = 'mean_diffs',
                null_iters = 1000,
                verbose = True,
                return_as_df = False,
                copy = False):
    return _sig_clusters(
        adata,
        obs_column_name,
        layer,
        filter_by_feature_groups,
        key_added,
        return_as_df,
        copy,
        sig_method = 'mean_diffs',
        null_iters = null_iters,
        verbose = verbose,
        pca_slot = None
    )










def _viper_similarity(adata,
                      nn=None,
                      ws=[4, 2],
                      alternative=['two-sided', 'greater', 'less'],
                      layer=None,
                      filter_by_feature_groups=None,
                      key_added='viper_similarity',
                      copy=False):
    if copy: adata = adata.copy()
    mat = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups).to_df()

    if np.min(mat.values.flatten())>=0:
        mat = pd.DataFrame(rankdata(mat,axis=1), index = mat.index, columns = mat.columns)
        row_sums = np.sum(~np.isnan(mat.values), axis=1).reshape(-1, 1)
        mat = pd.DataFrame(norm_ppf(mat/(row_sums+1)), index = mat.index, columns = mat.columns)

    mat[mat.isna()] =0 # will this work?

    xw = mat

    if nn == None:
        if alternative == 'greater':
            xw[xw < 0] = 0
        elif alternative == 'less' :
            xw[xw > 0] = 0

        if len(ws) == 1:
            xw = np.transpose(xw)/np.max(np.abs(mat), axis = 1)
            xw = np.sign(xw) * np.abs(xw) ** ws

        else:
            ws[1] = 1/(ws[1] - ws[0]) * np.log(1/0.9 -1)
            xw = np.sign(xw) *__sigT(np.abs(mat),ws[1],ws[0]) #why it's 1, 0 instead of 0,1

    else:
        if alternative == 'greater':
            xw = rankdata(-mat,axis=1)
            mat[xw > nn] = None
        elif alternative == 'less' :
            xw = rankdata(mat,axis=1)
            mat[xw > nn] = None
        else:
            xw = rankdata(mat,axis=1)
            mat[xw > nn/2 & xw <(len(xw) - nn/2 +1)] = None

    nes = np.sqrt(np.sum(xw**2, axis = 1))
    xw = xw.transpose()/np.sum(np.abs(xw),axis = 1)

    t2 = norm_ppf(rankdata(xw.transpose(), axis = 1)/(mat.shape[1]+1))
    vp = np.matmul(t2, xw)

    vp = vp * nes

    tmp = np.array([vp.values[np.triu_indices(vp.shape[0], 1)],vp.T.values[np.triu_indices(vp.shape[0], 1)]])
    tmp = np.sum(tmp * tmp ** 2, axis=0) / np.sum(tmp ** 2, axis=0)

    vp.values[np.triu_indices(vp.shape[0], 1)] = tmp
    vp = vp.T
    vp.values[np.triu_indices(vp.shape[0], 1)] = tmp
    vp.columns = vp.index

    adata.obsp[key_added] = vp

    if copy: return adata

def _aracne3_to_regulon(
    net_file,
    net_df,
    anno,
    MI_thres,
    regul_size,
    normalize_MI_per_regulon
):
    pd.options.mode.chained_assignment = None
    if net_df is None:
        net = pd.read_csv(net_file, sep='\t')
    else:
        net = net_df.copy()

    if anno is not None:
        if anno.shape[1] != 2 or not isinstance(anno, (pd.DataFrame, pd.Series, pd.Matrix)):
            raise ValueError("anno should contain two columns: 1-original symbol, 2-new symbol")

        anno.columns = ["old", "new"]

        ## Convert gene symbols:
        net['regulator.values'] = anno.set_index('old').loc[net['regulator.values'], 'new'].values
        net['target.values'] = anno.set_index('old').loc[net['target.values'], 'new'].values

    ## Network filtering
    net = net[net['mi.values'] > MI_thres]

    net.sort_values(by=['regulator.values','count.values','mi.values'],ascending=[True,False,False], inplace=True)
    net = net.groupby('regulator.values').head(regul_size)

    if normalize_MI_per_regulon:
        reg_max = net.groupby(['regulator.values'])['mi.values'].transform('max')
    else:
        reg_max = net['mi.values'].max()

    # print(reg_max)
    op = pd.DataFrame({
        'regulator': net['regulator.values'],
        'target' : net['target.values'],
        'mor': net['scc.values'],
        'likelihood': net['mi.values']/reg_max
        }
    )

    return op




def _adjust_p_values(p_values):
	# Helper function for *nes_to_pval* in module "tl"
    # correct p values with FDR
    _, corrected_p_values, _, _ = multitest.multipletests(p_values, method='fdr_bh')
    return corrected_p_values

def _nes_to_pval_df(dat_df, lower_tail=True, adjust = True, axs = 1, neg_log = False):
    """\
    Compute (adjusted) p-value associated to the viper-computed NES in a pd.DataFrame.

    Parameters
    ----------
    dat_df
        A pd.Series or pd.DataFrame containing protein activity (NES), pathways (NES) data or
        Stouffer-integrated NES data, where rows are observations/samples (e.g. cells or groups) and
        columns are features (e.g. proteins or pathways).
    lower_tail (default: True)
        If `True` (default), probabilities are P(X <= x)
        If `False`, probabilities are P(X > x)
    adjust (default: True)
        If `True`, returns adjusted p values using FDR Benjamini-Hochberg procedure.
        If `False`, does not adjust p values
    axs (default: 1)
        axis along which to perform the p-value correction (Used only if the input is a pd.DataFrame).
        Possible values are 0 or 1.

    Returns
    -------
    A pd.Series or pd.DataFrame objects of (adjusted) p-values.

    References
    ----------
    Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
        Journal of the Royal Statistical Society. Series B (Methodological), 57(1), 289–300.
        http://www.jstor.org/stable/2346101
    """

    if lower_tail == False:
        p_values_array = norm_sf(dat_df)

    elif lower_tail == True:
        p_values_array = 2 * norm_sf(np.abs(dat_df))

    else:
    	raise ValueError("'lower_tail' must be either True or False.")


    if dat_df.ndim == 1:
    # Calculate P values and corrected P values
        if adjust==True:
            _, p_values_array, _, _ = multitest.multipletests(p_values_array, method='fdr_bh')
            # Generate pd.DataFrames for (adjusted) p values
        p_values_df = pd.Series(p_values_array, index=dat_df.index)

    elif dat_df.ndim == 2:
    # Calculate P values and corrected P values
        if adjust==True:
            # correct p value
            p_values_array = np.apply_along_axis(_adjust_p_values, axis=axs, arr=p_values_array)

        # Generate pd.DataFrames for (adjusted) p values
        if isinstance(dat_df, pd.DataFrame):
            p_values_df = pd.DataFrame(p_values_array, index=dat_df.index, columns=dat_df.columns)
        elif isinstance(dat_df, np.ndarray):
            p_values_df = p_values_array
    else:
        raise ValueError("dat_df must have 1 or 2 dimensions.")

    if neg_log: p_values_df = -1 * np.log10(norm_sf(p_values_df))

    return p_values_df

def _nes_to_pval(
    adata,
    layer,
    key_added,
    lower_tail=True,
    adjust = True,
    axs = 1,
    neg_log = False,
    copy = False
):
    if copy: adata = adata.copy()

    if isinstance(adata, pd.DataFrame) or isinstance(adata, np.ndarray):
        adata[:] = _nes_to_pval_df(adata, lower_tail, adjust, axs, neg_log)
    elif(isinstance(adata, anndata.AnnData) or isinstance(adata, anndata._core.anndata.AnnData)):
        if layer is None:
            input_array = adata.X
        else:
            input_array = adata.layers[layer]

        transformed_array = _nes_to_pval_df(input_array, lower_tail, adjust, axs, neg_log)

        if key_added is not None:
            adata.layers[key_added] = transformed_array
        elif layer is not None:
            adata.layers[layer] = transformed_array
        else:
            adata.X = transformed_array
    else:
        raise Exception("adata must be anndata.AnnData, numpy.ndarray or pandas.DataFrame.")

    if copy: return adata















# from scipy.stats import rankdata
# from scipy.stats import spearmanr
# import time
#
# def Pearson_corr_single_X_col(x_diff, Y_diff, Sum_Y_diff_squared):
#     x_diff = x_diff[:, np.newaxis]
#     xY_diff = Y_diff * x_diff
#     Sum_xY_diff = np.sum(xY_diff, axis=0)
#     Sum_x_diff_squared = np.repeat(sum(x_diff**2), Y_diff.shape[1])
#     corr = Sum_xY_diff / np.sqrt(Sum_x_diff_squared * Sum_Y_diff_squared)
#     return corr
#
# def Pearson_corr_by_cols(X, Y):
#     X_mean = X.mean(axis=0, keepdims=True)
#     Y_mean = Y.mean(axis=0, keepdims=True)
#
#     X_diff = X - X_mean
#     Y_diff = Y - Y_mean
#
#     # For each Y_diff column we square and take the sum.
#     Sum_Y_diff_squared = np.sum(Y_diff**2, axis = 0)
#
#     # Make a correlation array for us to fill in
#     corr = np.apply_along_axis(Pearson_corr_single_X_col, 0, X_diff, Y_diff, Sum_Y_diff_squared)
#
#     return corr
#
# def Spearman_corr_by_cols(X, Y):
#     return Pearson_corr_by_cols(rankdata(X, axis=0), rankdata(Y, axis=0))
#
# def setdiff(X, Y):
#     return np.array([x for x in X if x not in Y])
#
# def get_spearman_network(norm_gex_data):
#     gex_data = norm_gex_data
#     s = time.time()
#     if isinstance(gex_data, anndata.AnnData):
#         gex_norm_df = gex_data.to_df()
#     else:
#         gex_norm_df = gex_data
#
#     TFs_CoTFs_array = np.array(pyviper.load.TFs() + pyviper.load.coTFs())
#     selected_columns = np.intersect1d(TFs_CoTFs_array, gex_norm_df.columns.values)
#
#     spearman_corr_array = Spearman_corr_by_cols(
#         gex_norm_df[selected_columns].values,
#         gex_norm_df.drop(columns=selected_columns).values
#     )
#
#     spearman_net_table = pd.DataFrame(spearman_corr_array)
#     spearman_net_table.columns = selected_columns
#     spearman_net_table.index = setdiff(gex_norm_df.columns.values, selected_columns)
#     spearman_net_table = spearman_net_table.reset_index()
#     spearman_net_table.rename(columns={'index': 'target'}, inplace=True)
#     spearman_net_table = pd.melt(spearman_net_table, id_vars='target', var_name="regulator", value_name='mor')
#     spearman_net_table.dropna(subset=['mor'], inplace=True)
#     spearman_net_table['likelihood'] = np.abs(spearman_net_table['mor'].values)
#     spearman_net_table = spearman_net_table[['regulator', 'target', 'mor', 'likelihood']]
#     spearman_net = pyviper.Interactome("spearman_net", spearman_net_table)
#     e = time.time()
#     print (e-s)
#
#     return spearman_net
