### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import rankdata
from scipy.stats import norm
from statsmodels.stats import multitest
from ._filtering_funcs import *
from ._filtering_funcs import _get_anndata_filtered_by_feature_group
from ._helpers import _adjust_p_values

### ---------- EXPORT LIST ----------
__all__ = []

# ------------------------ SCANPY TOOLS PYTHER WRAPPERS -----------------------
def pca(adata,
        *,
        layer=None,
        filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
        **kwargs):
    """\
    A wrapper for the scanpy function sc.tl.pca.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    layer (default: None)
        The layer to use as input data.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    **kwargs
        Arguments to provide to the sc.tl.pca function.
    """
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    print(adata_filt.shape)

    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]

def dendrogram(adata,
               *,
               groupby,
               key_added=None,
               layer=None,
               filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
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
    filter_by_feature_groups (default: None)
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

def stouffer(adata, obs_column_name, layer = None, filter_by_feature_groups=None):
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
    filter_by_feature_groups (default: None)
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
        A pd.Series or pd.DataFrame containing protein activity (NES), pathways (NES) data or
        Stouffer-integrated NES data, where rows are observations/samples (e.g. cells or groups) and
        columns are features (e.g. proteins or pathways).
    adjust (default: True)
        If `True`, returns adjusted p values using FDR Benjamini-Hochberg procedure.
        If `False`, does not adjust p values

    Returns
    -------
    A pd.Series or pd.DataFrame objects of (adjusted) p-values.

    References
    ----------
    Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
        Journal of the Royal Statistical Society. Series B (Methodological), 57(1), 289–300.
        http://www.jstor.org/stable/2346101
    """

    p_values_array = 2 * norm.sf(np.abs(dat_df))

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
            p_values_array = np.apply_along_axis(_adjust_p_values, axis=1, arr=p_values_array)

        # Generate pd.DataFrames for (adjusted) p values
        p_values_df = pd.DataFrame(p_values_array, index=dat_df.index, columns=dat_df.columns)
    else:
        raise ValueError("dat_df must have 1 or 2 dimensions.")

    return p_values_df

def _generate_interactome_from_pax_data(pax_data,
                                        interactome_name="vpmat",
                                        n_top=50,
                                        is_symmetric=True):

    if isinstance(pax_data, anndata.AnnData):
        vpmat = pax_data.to_df()
    elif isinstance(pax_data, pd.DataFrame):
        vpmat = pax_data
    else:
        raise ValueError("pax_data must be anndata.AnnData or pd.DataFrame.")


    n_cells = vpmat.shape[0]
    n_mrs = vpmat.shape[1]
    cell_names = vpmat.index

    # For each sample, we calculate index arragement that would sort the vector
    sorted_order_array = np.argsort(-vpmat.values)
    # We then get the MRs ranked for each sample by indexing with this sorted order
    mrs_ranked_array = np.array(vpmat.columns)[sorted_order_array]

    if is_symmetric:
        # Slice the top n_top/2 rows and bottom n_top/2 rows
        # Get the top 25 and bottom 25 rows
        selected_column_indices = list(range(0,n_top_half)) + list(range(n_mrs-n_top_half,n_mrs))
        cell_i_mor = np.concatenate((np.ones(n_top_half), np.full(n_top_half, -1)))
    else:
        # Slice the top n_top rows
        selected_column_indices = list(range(0,n_top))
        cell_i_mor = np.ones(n_top)

    top_mrs_ranked_array = mrs_ranked_array[:, selected_column_indices]
    top_mrs_ranked_array_1D = top_mrs_ranked_array.flatten()

    regulator = np.repeat(cell_names, n_top)
    target = top_mrs_ranked_array_1D
    mor = np.tile(cell_i_mor, n_cells)

    net_table = pd.DataFrame({
            'regulator': regulator,
            'target': top_mrs_ranked_array_1D,
            'mor': mor,
            'likelihood': 1
        })

    return pyther.Interactome(interactome_name, net_table)



def OncoMatch(pax_data_to_test,
              pax_data_for_cMRs,
              tcm_size = 50,
              both_ways = False,
              om_max_NES_threshold = 30,
              om_min_logp_threshold = 0,
              enrichment = 'area',
              key_added = 'om'):
    """\
    The OncoMatch function that computes -log10 p-values for each sample in pax_data_to_test
    of the MRs of each sample in pax_data_for_cMRs.

    Parameters
    ----------
    pax_data_to_test
        An anndata.AnnData or pd.DataFrame containing protein activity (NES),
        where rows are observations/samples (e.g. cells or groups) and
        columns are features (e.g. proteins or pathways).
    pax_data_for_cMRs
        An anndata.AnnData or pd.DataFrame containing protein activity (NES),
        where rows are observations/samples (e.g. cells or groups) and
        columns are features (e.g. proteins or pathways).
    tcm_size (default: 50)
        Number of top MRs from each sample to use to compute regulators.
    both_ways (default: False)
        Whether to also use the candidate MRs of pax_data_to_test to compute
        NES for the samples in pax_data_for_cMRs, and then average.
    om_max_NES_threshold (default: 30)
        The maximum NES scores before using a cutoff.
    om_min_logp_threshold (default: 0)
        The minimum logp value threshold, such that all logp values smaller than
        this value are set to 0.
    enrichment (default: 'area')
        The method of compute enrichment. 'area' or 'narnea'
    key_added (default: 'om')
        The slot in pax_data_to_test.obsm to store the oncomatch results.
    Returns
    -------
    A pd.DataFrame objects of -log10 p-values with shape n_samples in
    pax_data_to_test by n_samples pax_data_for_cMRs.

    References
    ----------
    Alvarez, M. J. et al. A precision oncology approach to the pharmacological
    targeting of mechanistic dependencies in neuroendocrine tumors. Nat Genet 50,
    979–989, doi:10.1038/s41588-018-0138-4 (2018).

    Alvarez, M. J. et al. Reply to ’H-STS, L-STS and KRJ-I are not authentic GEPNET
    cell lines’. Nat Genet 51, 1427–1428, doi:10.1038/s41588-019-0509-5 (2019).
    """

    if isinstance(pax_data_to_test, anndata.AnnData):
        vpmat_to_test = pax_data_to_test.to_df()
    elif isinstance(pax_data_to_test, pd.DataFrame):
        vpmat_to_test = pax_data_to_test
    else:
        raise ValueError("pax_data_to_test must be anndata.AnnData or pd.DataFrame.")

    if isinstance(pax_data_for_cMRs, anndata.AnnData):
        vpmat_for_cMRs = pax_data_for_cMRs.to_df()
    elif isinstance(pax_data_for_cMRs, pd.DataFrame):
        vpmat_for_cMRs = pax_data_for_cMRs
    else:
        raise ValueError("pax_data_for_cMRs must be anndata.AnnData or pd.DataFrame.")

    # Compute intersection of regulons
    regs_in_common = np.intersect1d(vpmat_to_test.columns.values, vpmat_for_cMRs.columns.values)
    vpmat_to_test = vpmat_to_test[regs_in_common]
    vpmat_for_cMRs = vpmat_for_cMRs[regs_in_common]


    if both_ways is True:
        interactome_test = _generate_interactome_from_pax_data(
            vpmat_to_test,
            interactome_name = "vpmat_to_test",
            n_top = tcm_size
        )
        interactome_cMRs = _generate_interactome_from_pax_data(
            vpmat_for_cMRs,
            interactome_name = "vpmat_to_test",
            n_top = tcm_size
        )

        om_t = pyther.viper(gex_data = anndata.AnnData(vpmat_to_test, dtype='float64'),
                            interactome = interactome_cMRs,
                            enrichment = enrichment,
                            min_targets = 0,
                            output_as_anndata=False,
                            verbose = False)
        if enrichment == 'narnea': om_t = om_t['nes']

        om_q = pyther.viper(gex_data = anndata.AnnData(vpmat_for_cMRs, dtype='float64'),
                            interactome = interactome_test,
                            enrichment = enrichment,
                            min_targets = 0,
                            output_as_anndata=False,
                            verbose = False)
        if enrichment == 'narnea': om_q = om_q['nes']

        # Replace NaN (missing) values with 0 in om_t
        om_t[np.isnan(om_t)] = 0

        # Replace NaN (missing) values with 0 in om_q
        om_q[np.isnan(om_q)] = 0

        # Clip values greater than om_max_NES_threshold in om_t
        om_t = np.where(om_t > om_max_NES_threshold, om_max_NES_threshold, om_t)

        # Clip values greater than om_max_NES_threshold in om_q
        om_q = np.where(om_q > om_max_NES_threshold, om_max_NES_threshold, om_q)

        # Tranpose om_q so it has the same shape as om_t
        om_q = np.transpose(om_q)

        # Average them
        om = (om_t+om_q)/2
    else:
        interactome_cMRs = _generate_interactome_from_pax_data(
            vpmat_for_cMRs,
            interactome_name = "vpmat_to_test",
            n_top = tcm_size
        )
        om = pyther.viper(gex_data = anndata.AnnData(vpmat_to_test, dtype='float64'),
                          interactome = interactome_cMRs,
                          enrichment = enrichment,
                          min_targets = 0,
                          output_as_anndata=False,
                          verbose = False)
        if enrichment == 'narnea': om = om['nes']

        # Replace NaN (missing) values with 0 in om
        om[np.isnan(om)] = 0

        # Clip values greater than om_max_NES_threshold in om
        om = np.where(om > om_max_NES_threshold, om_max_NES_threshold, om)

    # Compute p_values
    # cond = om < 7
    # om[cond] = 1 - norm.cdf(om[cond]) # accurate (i.e. same as R) when NES scores are small (returns 0.0 when given big NES, e.g. 10)
    # om[~cond] = norm.logcdf(om[~cond])*-1 # accurate (i.e. same as R) when NES scores are big (e.g. 10)

    om = pd.DataFrame(om, index = vpmat_to_test.index, columns = vpmat_for_cMRs.index)
    om = nes_to_pval_df(om)

    # Log transform
    om = -np.log10(om)

    # Clip values smaller than om_min_logp_threshold in om
    om = np.where(om < om_min_logp_threshold, 0, om)

    om = pd.DataFrame(om, index = vpmat_to_test.index, columns = vpmat_for_cMRs.index)

    pax_data_to_test.obsm[key_added] = om
