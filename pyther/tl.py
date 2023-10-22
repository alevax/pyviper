### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import norm
from ._filtering_funcs import *
from ._filtering_funcs import _get_anndata_filtered_by_feature_group
from ._helpers import _adjust_p_values

### ---------- EXPORT LIST ----------
__all__ = []

# ------------------------ SCANPY TOOLS PYTHER WRAPPERS -----------------------
def pca(adata,
        *,
        layer=None,
        filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
        **kwargs):
    """\
    A wrapper for the scanpy function sc.tl.pca.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    layer (default: None)
        The layer to use as input data.
    filter_by_feature_groups (default: ["tfs", "cotfs"])
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    **kwargs
        Arguments to provide to the sc.tl.pca function.
    """
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)

    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]

def dendrogram(adata,
               *,
               groupby,
               key_added=None,
               layer=None,
               filter_by_feature_groups=["tfs", "cotfs"], # ["tfs", "cotfs", "sig", "surf"],
               **kwargs):
    """\
    A wrapper for the scanpy function sc.tl.dendrogram.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    key_added (default: None)
        The key in adata.uns where the dendrogram should be stored.
    layer (default: None)
        The layer to use as input data.
    filter_by_feature_groups (default: ["tfs", "cotfs"])
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    **kwargs
        Arguments to provide to the sc.tl.dendrogram function.
    """
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    if key_added is None:
        key_added = f'dendrogram_{groupby}'
    sc.tl.dendrogram(adata_filt, groupby, **kwargs, key_added = key_added)
    adata.uns[key_added] = adata_filt.uns[key_added]

def stouffer(adata, obs_column_name, layer = None, filter_by_feature_groups=["tfs", "cotfs"]):
    """\
    Compute a stouffer signature on each of your clusters in an anndata object.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    obs_column_name
        The name of the column of observations to use as clusters.
    layer (default: None)
        The layer to use as input data to compute the signatures.
    filter_by_feature_groups (default: ["tfs", "cotfs"])
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".

    Returns
    -------
    A new anndata object containing cluster stouffer signatures.
    """
    adata = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    if layer is None:
        dat_df = pd.DataFrame(adata.X,
                              index=adata.obs_names,
                              columns=adata.var_names)
    else:
        dat_df = pd.DataFrame(adata.layers[layer],
                              index=adata.obs_names,
                              columns=adata.var_names)
    cluster_vector = adata.obs[obs_column_name]
    result_df = stouffer_clusters_df(dat_df, cluster_vector)
    adata.uns['stouffer'] = result_df
    return adata
    # return mat_to_anndata(result_df)

def stouffer_clusters_df(dat_df, cluster_vector):
    """\
    Compute a stouffer signature on each of your clusters from a DataFrame.

    Parameters
    ----------
    dat_df
        A pandas dataframe containing input data.
    cluster_vector
        A cluster vector corresponding to observations in the pd.DataFrame.

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



def nes_to_pval_df(dat_df,adjust=True):
    """\
    Compute (adjusted) p-value associated to the viper-computed NES in a pd.DataFrame.

    Parameters
    ----------
    dat_df
        A pandas dataframe containing protein activity (NES), pathways (NES) data or 
        Stouffer-integrated NES data, where rows are observations/samples (e.g. cells or groups) and 
        columns are features (e.g. proteins or pathways).
    adjust
        If `True`, returns adjusted p values using FDR Benjamini-Hochberg procedure.
        If `False`, does not adjust p values
    
    Returns
    -------
    A pd.DataFrame objects of (adjusted) p-values.

    References
    ----------
    Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. 
        Journal of the Royal Statistical Society. Series B (Methodological), 57(1), 289–300. 
        http://www.jstor.org/stable/2346101
    """
    
    p_values_array = 2 * norm.sf(np.abs(dat_df))

    if dat_df.ndim == 1:
    # Calculate P values and corrected P values
 
    elif dat_df.ndim == 2:
    # Calculate P values and corrected P values
        if adjust==True:
            # correct p value 
            p_values_array = np.apply_along_axis(_adjust_p_values, axis=1, arr=dat_df, adjust=True)
        elif adjust==False:
            # do not correct for p value
            p_values_array = np.apply_along_axis(_adjust_p_values, axis=1, arr=dat_df, adjust=False)
        else:
            raise ValueError("Parameter adjust must be either True or False [bool]") 
    else:
        raise ValueError("dat_df must have 1 or 2 dimensions.") 

    # Generate pd.DataFrames for p values and corrected p values
    p_values_df = pd.DataFrame(p_values_array, index=dat_df.index, columns=dat_df.columns)

    return p_values_df 


    