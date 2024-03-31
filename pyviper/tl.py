### ---------- IMPORT DEPENDENCIES ----------
from ._tl import _pca, _dendrogram, _oncomatch, _find_top_mrs, _path_enr

### ---------- EXPORT LIST ----------
__all__ = []

# ------------------------ SCANPY TOOLS PYVIPER WRAPPERS -----------------------
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
    _pca(adata,
         layer,
         filter_by_feature_groups,
         **kwargs)

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
    _dendrogram(adata,
                groupby,
                key_added,
                layer,
                filter_by_feature_groups,
                **kwargs)

def oncomatch(pax_data_to_test,
               pax_data_for_cMRs,
               tcm_size = 50,
               both_ways = False,
               lower_tail = True,
               om_max_NES_threshold = 30,
               om_min_logp_threshold = 0,
               enrichment = 'aREA',
               key_added = 'om',
               return_as_df = False,
               copy = False):
    """\
    The OncoMatch algorithm[1] assesses the overlap in differentially active MR
    proteins between two sets of samples (e.g. to validate GEMMs as effective
    models of human tumor samples). It does so by computing -log10 p-values for
    each sample in pax_data_to_test of the MRs of each sample in pax_data_for_cMRs.

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
    lower_tail: default (True)
    	If `True` (default), probabilities are P(X <= x)
    	If `False`, probabilities are P(X > x)
    om_max_NES_threshold (default: 30)
        The maximum NES scores before using a cutoff.
    om_min_logp_threshold (default: 0)
        The minimum logp value threshold, such that all logp values smaller than
        this value are set to 0.
    enrichment (default: 'aREA')
        The method of compute enrichment. 'aREA' or 'NaRnEA'
    key_added (default: 'om')
        The slot in pax_data_to_test.obsm to store the oncomatch results.
    return_as_df (default: False)
        Instead of adding the OncoMatch DataFrame to pax_data_to_test.obsm,
        return it directly.
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    When copy is False, stores a pd.DataFrame objects of -log10 p-values with
    shape (n_samples in pax_data_to_test, n_samples in pax_data_for_cMRs) in
    pax_data_to_test.obsm[key_added]. When copy is True, a copy of the AnnData
    is returned with these pd.DataFrames stored. When return_as_df is True,
    the OncoMatch DataFrame alone is directly returned by the function.

    References
    ----------
    Alvarez, M. J. et al. A precision oncology approach to the pharmacological
    targeting of mechanistic dependencies in neuroendocrine tumors. Nat Genet 50,
    979–989, doi:10.1038/s41588-018-0138-4 (2018).

    Alvarez, M. J. et al. Reply to ’H-STS, L-STS and KRJ-I are not authentic GEPNET
    cell lines’. Nat Genet 51, 1427–1428, doi:10.1038/s41588-019-0509-5 (2019).
    """
    return _oncomatch(pax_data_to_test,
                      pax_data_for_cMRs,
                      tcm_size,
                      both_ways,
                      lower_tail,
                      om_max_NES_threshold,
                      om_min_logp_threshold,
                      enrichment,
                      key_added,
                      return_as_df,
                      copy)

def find_top_mrs(adata,
                 pca_slot = "X_pca",
                 obs_column_name = None,
                 layer = None,
                 N = 50,
                 both = True,
                 method = "stouffer",
                 key_added = "mr",
                 filter_by_feature_groups=None,
                 rank=False,
                 filter_by_top_mrs = False,
                 return_as_df = False,
                 copy = False):
    """\
    Identify the top N master regulator proteins in a VIPER AnnData object

    Parameters
    ----------
    adata
        An anndata object containing a distance object in adata.obsp.
    pca_slot
        The slot in adata.obsm where a PCA is stored. Only required when method
        is "spearman".
    obs_column_name
        The name of the column of observations in adata to use as clusters, or a
        cluster vector corresponding to observations. Required when method is
        "mwu" or "spearman".
    N (default: 50)
        The number of MRs to return
    both (default: True)
        Whether to return both the top N and bottom N MRs (True) or just the
        top N (False).
    method (default: "stouffer")
        The method used to compute a signature to identify the top candidate
        master regulators (MRs). The options come from functions in pyviper.pp.
        Choose between "stouffer", "mwu", or "spearman".
    key_added (default: "mr")
        The name of the slot in the adata.var to store the output.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    rank (default: False)
        When False, a column is added to var with identified MRs labeled as
        "True", while all other proteins are labeled as "False". When True, top
        MRs are labeled N,N-1,N-2,...,1, bottom MRs are labeled -N,-N-1,-N-2,
        ...,-1, and all other proteins are labeled 0. Higher rank means greater
        activity, while lower rank means less.
    filter_by_top_mrs (default: False)
        Whether to filter var to only the top MRs in adata
    return_as_df (default: False)
        Returns a pd.DataFrame of the top MRs per cluster
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    Add a column to adata.var[key_added] or, when clusters given, adds multiple
    columns (e.g. key_added_clust1name, key_added_clust2name, etc) to adata.var.
    If copy, returns a new adata transformed by this function. If return_as_df,
    returns a DataFrame.
    """
    # Feature where you can choose a method, e.g. MWU Test instead of Stouffer signature
    # scipy.stats.mannwhitneyu
    return _find_top_mrs(
        adata,
        pca_slot,
        obs_column_name,
        layer,
        N,
        both,
        method,
        key_added,
        filter_by_feature_groups,
        rank,
        filter_by_top_mrs,
        return_as_df,
        copy
    )

def path_enr(gex_data,
             pathway_interactome,
             layer=None,
             eset_filter=True,
             method=None,  # [None, "scale", "rank", "mad", "ttest"],
             enrichment='aREA',  # [None, 'area','narnea'],
             mvws=1,
             njobs=1,
             batch_size=10000,
             verbose=True,
             output_as_anndata=True,
             transfer_obs=True,
             store_input_data=True
             ):
    """\
    Run the variation of VIPER that is specific to pathway enrichment analysis:
    a single interactome and min_targets is set to 0.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy).
    pathway_interactome
        An object of class Interactome or one of the following strings that
        corresponds to msigdb regulons: "c2", "c5", "c6", "c7", "h".
    layer (default: None)
        The layer in the anndata object to use as the gene expression input.
    eset_filter (default: False)
        Whether to filter out genes not present in the interactome (True) or to
        keep this biological context (False). This will affect gene rankings.
    method (default: None)
        A method used to create a gene expression signature from gex_data.X. The
        default of None is used when gex_data.X is already a gene expression
        signature. Alternative inputs include "scale", "rank", "doublerank",
        "mad", and "ttest".
    enrichment (default: 'aREA')
        The algorithm to use to calculate the enrichment. Choose betweeen
        Analytical Ranked Enrichment Analysis (aREA) and Nonparametric
        Analytical Rank-based Enrichment Analysis (NaRnEA) function. Default ='aREA',
        alternative = 'NaRnEA'.
    mvws (default: 1)
        (A) Number indicating either the exponent score for the metaViper weights.
        These are only applicable when enrichment = 'aREA' and are not used when
        enrichment = 'NaRnEA'. Roughly, a lower number (e.g. 1) results in
        networks being treated as a consensus network (useful for multiple
        networks of the same celltype with the same epigenetics), while a higher
        number (e.g. 10) results in networks being treated as separate (useful
        for multiple networks of different celltypes with different epigenetics).
        (B) The name of a column in gex_data that contains the manual assignments
        of samples to networks using list position or network names.
        (C) "auto": assign samples to networks based on how well each
        network allows for sample enrichment.
    njobs (default: 1)
        Number of cores to distribute sample batches into.
    batch_size (default: 10000)
        Maximum number of samples to process at once. Set to None to split all
        samples across provided `njobs`.
    verbose (default: True)
        Whether extended output about the progress of the algorithm should be
        given.
    output_as_anndata (default: True)
        Way of delivering output.
    transfer_obs (default: True)
        Whether to transfer the observation metadata from the input anndata to
        the output anndata. Thus, not applicable when output_as_anndata==False.
    store_input_data (default: True)
        Whether to store the input anndata in an unstructured data slot (.uns) of
        the output anndata. Thus, not applicable when output_as_anndata==False.
        If input anndata already contains 'gex_data' in .uns, the input will
        assumed to be protein activity and will be stored in .uns as 'pax_data'.
        Otherwise, the data will be stored as 'gex_data' in .uns.

    Returns
    -------
    Returns an AnnData object containing the pathways. When store_input_data,
    the input gex_data AnnData is stored within the dataframe.
    """
    return _path_enr(gex_data,
                     pathway_interactome,
                     layer,
                     eset_filter,
                     method,
                     enrichment,
                     mvws,
                     njobs,
                     batch_size,
                     verbose,
                     output_as_anndata,
                     transfer_obs,
                     store_input_data)
