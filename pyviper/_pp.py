import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from scipy.stats import rankdata, norm
from statsmodels.stats import multitest
from ._filtering_funcs import _get_anndata_filtered_by_feature_group
from tqdm.auto import tqdm
import os

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
    key_added = None
):
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


def __sigT(x, slope = 20, inflection = 0.5):
    return (1 - 1/(1 + np.exp(slope * (x - inflection))))

def _stouffer_clusters_df(dat_df, cluster_vector):
    """\
    Compute a stouffer signature on each of your clusters from a DataFrame.

    Parameters
    ----------
    dat_df

    cluster_vector


    Returns
    -------
    A new pd.DataFrame containing cluster stouffer signatures.
    """
    # Ensure cluster_vector has the same number of samples as rows in dat_df
    if len(cluster_vector) != dat_df.shape[0]:
        raise ValueError("Cluster vector length does not match the number of rows in the DataFrame.")

    # Convert the DataFrame to a NumPy array
    dat_array = dat_df.to_numpy()

    # Find unique clusters and initialize arrays to store Stouffer scores
    unique_clusters, cluster_indices = np.unique(cluster_vector, return_inverse=True)
    n_clusters = len(unique_clusters)
    n_genes = dat_df.shape[1]
    stouffer_scores = np.zeros((n_clusters, n_genes))

    # Calculate the denominator for Stouffer scores for each cluster
    cluster_sizes = np.bincount(cluster_indices)
    sqrt_cluster_sizes = np.sqrt(cluster_sizes)

    # Calculate Stouffer scores for each cluster and gene
    for i in range(n_clusters):
        cluster_mask = (cluster_indices == i)
        cluster_data = dat_array[cluster_mask]
        stouffer_scores[i, :] = np.sum(cluster_data, axis=0) / sqrt_cluster_sizes[i]

    # Create a DataFrame from the computed Stouffer scores
    result_df = pd.DataFrame(stouffer_scores, index=unique_clusters, columns=dat_df.columns)

    return result_df

def _stouffer_clusters_adata(adata,
                             obs_column_name = None,
                             layer = None,
                             filter_by_feature_groups = None):
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

    result_df = _stouffer_clusters_df(dat_df, cluster_vector)
    return result_df

def _stouffer(adata,
              obs_column_name = None,
              layer = None,
              filter_by_feature_groups = None,
              key_added = 'stouffer'):
    result_df = _stouffer_clusters_adata(adata,
                                         obs_column_name,
                                         layer,
                                         filter_by_feature_groups).T
    if obs_column_name is None:
        result_df.columns = [key_added]
    else:
        result_df.columns = key_added + "_" + result_df.columns
    adata.var = pd.concat([adata.var, result_df], axis=1, join='inner')

def _viper_similarity(adata,
                         nn = None,
                         ws = [4, 2],
                         alternative=['two-sided','greater','less'],
                         layer=None,
                         filter_by_feature_groups=None,
                         key_added = 'viper_similarity'):
    mat = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups).to_df()

    if np.min(mat)>=0 :
        mat = rankdata(mat,axis=1)
        mat = norm.ppf(mat/(np.sum(mat.isna()==False,axis = 1)+1))

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

    t2 = norm.ppf(rankdata(xw.transpose(), axis = 1)/(mat.shape[1]+1))
    vp = np.matmul(t2, xw)

    vp = vp * nes

    tmp = np.array([vp.values[np.triu_indices(vp.shape[0], 1)],vp.T.values[np.triu_indices(vp.shape[0], 1)]])
    tmp = np.sum(tmp * tmp ** 2, axis=0) / np.sum(tmp ** 2, axis=0)

    vp.values[np.triu_indices(vp.shape[0], 1)] = tmp
    vp = vp.T
    vp.values[np.triu_indices(vp.shape[0], 1)] = tmp
    vp.columns = vp.index

    adata.obsp[key_added] = vp

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

# def _nes_to_neg_log(adata, layer = None, key_added = None):
#     if isinstance(adata, pd.DataFrame) or isinstance(adata, np.ndarray):
#         adata[:] = -1 * np.log10(norm.sf(adata))
#     elif(isinstance(adata, anndata.AnnData) or isinstance(adata, anndata._core.anndata.AnnData)):
#         if layer is None:
#             input_array = adata.X
#         else:
#             input_array = adata.layers[layer]
#
#         transformed_array = -1*np.log10(norm.sf(input_array))
#
#         if key_added is not None:
#             adata.layers[key_added] = transformed_array
#         elif layer is not None:
#             adata.layers[layer] = transformed_array
#         else:
#             adata.X = transformed_array
#     else:
#         raise Exception("adata must be anndata.AnnData, numpy.ndarray or pandas.DataFrame.")
#






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
        p_values_array = norm.sf(dat_df)

    elif lower_tail == True:
        p_values_array = 2 * norm.sf(np.abs(dat_df))

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

    if neg_log: p_values_df = -1 * np.log10(norm.sf(p_values_df))

    return p_values_df

def _nes_to_pval(adata, layer, key_added, lower_tail=True, adjust = True, axs = 1, neg_log = False):
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
