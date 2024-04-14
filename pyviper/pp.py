### ---------- IMPORT DEPENDENCIES ----------
from ._pp import _rank_norm, _stouffer, _mwu, _spearman, _viper_similarity, _aracne3_to_regulon, _nes_to_pval, _mad_from_R, _median
from ._corr_distance import corr_distance
from ._rep_subsample_funcs import _representative_subsample_anndata
from ._metacell_funcs import _representative_metacells_multiclusters
from ._translate import translate

### ---------- EXPORT LIST ----------
__all__ = []

def rank_norm(
    adata,
    NUM_FUN=_median,
    DEM_FUN = _mad_from_R,
    layer = None,
    key_added = None,
    copy = False
):
    """\
    Compute a double rank normalization on an anndata, np.array, or pd.DataFrame.

    Parameters
    ----------
    adata
        Data stored in an anndata object, np.array or pd.DataFrame.
    NUM_FUN (default: np.median)
        The first function to be applied across each column.
    DEM_FUN (default: _mad_from_R)
        The second function to be applied across each column.
    layer (default: None)
        For an anndata input, the layer to use. When None, the input layer is
        anndata.X.
    key_added (default: None)
        For an anndata input, the name of the layer where to store. When None,
        this is anndata.X.
    copy (default: False)
        Whether to return a rank-transformed copy (True) or to instead transform
        the original input (False).

    Returns
    -------
    When copy = False, saves the input data as a double rank transformed version.
    When copy = True, return a double rank transformed version of the input data.
    """
    return _rank_norm(
        adata,
        NUM_FUN,
        DEM_FUN,
        layer,
        key_added,
        copy
    )

def stouffer(adata,
             obs_column_name = None,
             layer = None,
             filter_by_feature_groups = None,
             key_added = 'stouffer',
             compute_pvals = True,
             null_iters = 1000,
             verbose = True,
             return_as_df = False,
             copy = False):
    """\
    Compute a stouffer signature on each of your clusters in an anndata object.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata
        object, or a pandas dataframe containing input data.
    obs_column_name
        The name of the column of observations in adata to use as clusters, or a
        cluster vector corresponding to observations.
    layer (default: None)
        The layer to use as input data to compute the signatures.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    key_added (default: 'stouffer')
        The slot in adata.uns to store the stouffer signatures.
    compute_pvals (default: True)
        Whether to compute a p-value for each score to return in the results.
    null_iters (default: 1000)
        The number of iterations to use to compute a null model to assess the
        p-values of each of the stouffer scores.
    verbose (default: True)
        Whether to provide additional output during the execution of the function.
    return_as_df (default: False)
        If True, returns the stouffer signature in a pd.DataFrame. If False,
        stores it in adata.var[key_added].
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    When return_as_df is False, adds the cluster stouffer signatures to
    adata.var[key_added]. When return_as_df is True, returns as pd.DataFrame.
    """
    return _stouffer(adata,
                     obs_column_name,
                     layer,
                     filter_by_feature_groups,
                     key_added,
                     compute_pvals,
                     null_iters,
                     verbose,
                     return_as_df,
                     copy)

def mwu(adata,
        obs_column_name = None,
        layer = None,
        filter_by_feature_groups = None,
        key_added = 'mwu',
        compute_pvals = True,
        verbose = True,
        return_as_df = False,
        copy = False):
    """\
    Compute a Mann-Whitney U-Test signature on each of your clusters in an
    anndata object.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata
        object, or a pandas dataframe containing input data.
    obs_column_name
        The name of the column of observations in adata to use as clusters, or a
        cluster vector corresponding to observations.
    layer (default: None)
        The layer to use as input data to compute the signatures.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    key_added (default: 'mwu')
        The slot in adata.uns to store the MWU signatures.
    compute_pvals (default: True)
        Whether to compute a p-value for each score to return in the results.
    verbose (default: True)
        Whether to provide additional output during the execution of the function.
    return_as_df (default: False)
        If True, returns the MWU signature in a pd.DataFrame. If False,
        stores it in adata.var[key_added].
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    When return_as_df is False, adds the cluster MWU signatures to
    adata.var[key_added]. When return_as_df is True, returns as pd.DataFrame.
    """
    return _mwu(adata,
                obs_column_name,
                layer,
                filter_by_feature_groups,
                key_added,
                compute_pvals,
                verbose,
                return_as_df,
                copy)

def spearman(adata,
             pca_slot = "X_pca",
             obs_column_name = None,
             layer = None,
             filter_by_feature_groups = None,
             key_added = 'stouffer',
             compute_pvals = True,
             null_iters = 1000,
             verbose = True,
             return_as_df = False,
             copy = False):
    """\
    Compute spearman correlation between each gene product and the cluster
    centroids along with the statistical significance for each of your clusters
    in an anndata object.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata
        object, or a pandas dataframe containing input data.
    pca_slot
        The slot in adata.obsm where a PCA is stored.
    obs_column_name
        The name of the column of observations in adata to use as clusters, or a
        cluster vector corresponding to observations.
    layer (default: None)
        The layer to use as input data to compute the correlation.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    key_added (default: 'spearman')
        The slot in adata.uns to store the spearman correlation.
    compute_pvals (default: True)
        Whether to compute a p-value for each score to return in the results.
    null_iters (default: 1000)
        The number of iterations to use to compute a null model to assess the
        p-values of each of the spearman scores.
    verbose (default: True)
        Whether to provide additional output during the execution of the function.
    return_as_df (default: False)
        If True, returns the spearman signature in a pd.DataFrame. If False,
        stores it in adata.var[key_added].
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    When return_as_df is False, adds the cluster spearman correlation to
    adata.var[key_added]. When return_as_df is True, returns as pd.DataFrame.
    """
    return _spearman(adata,
                     pca_slot,
                     obs_column_name,
                     layer,
                     filter_by_feature_groups,
                     key_added,
                     compute_pvals,
                     null_iters,
                     verbose,
                     return_as_df,
                     copy)


def viper_similarity(adata,
                     nn = None,
                     ws = [4, 2],
                     alternative=['two-sided','greater','less'],
                     layer=None,
                     filter_by_feature_groups=None,
                     key_added = 'viper_similarity',
                     copy = False):
    """\
    Compute the similarity between the columns of a VIPER-predicted activity or
    gene expression matrix. While following the same concept as the two-tail
    Gene Set Enrichment Analysis (GSEA)[1], it is based on the aREA algorithm[2].

    If ws is a single number, weighting is performed using an exponential function.
    If ws is a 2 numbers vector, weighting is performed with a symmetric sigmoid
    function using the first element as inflection point and the second as trend.

    Parameters
    ----------
    adata
        An anndata.AnnData containing protein activity (NES), where rows are
        observations/samples (e.g. cells or groups) and columns are features
        (e.g. proteins or pathways).
    nn (default: None)
        Optional number of top regulators to consider for computing the similarity
    ws (default: [4, 2])
        Number indicating the weighting exponent for the signature, or vector of
        2 numbers indicating the inflection point and the value corresponding to
        a weighting score of .1 for a sigmoid transformation, only used if nn is
        ommited.
    alternative (default: 'two-sided')
        Character string indicating whether the most active (greater), less
        active (less) or both tails (two.sided) of the signature should be used
        for computing the similarity.
    layer (default: None)
        The layer to use as input data to compute the signatures.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".
    key_added (default: "viper_similarity")
        The name of the slot in the adata.obsp to store the output.
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    Saves a signature-based distance numpy.ndarray in adata.obsp[key_added].

    References
    ----------
    [1] Julio, M. K. -d. et al. Regulation of extra-embryonic endoderm stem cell
    differentiation by Nodal and Cripto signaling. Development 138, 3885-3895 (2011).
    [2] Alvarez, M. J., Shen, Y., Giorgi, F. M., Lachmann, A., Ding, B. B., Ye, B. H.,
    & Califano, A. (2016). Functional characterization of somatic mutations in
    cancer using network-based inference of protein activity. Nature genetics,
    48(8), 838-847.
    """
    return _viper_similarity(
        adata,
        nn,
        ws,
        alternative,
        layer,
        filter_by_feature_groups,
        key_added,
        copy
    )

def aracne3_to_regulon(
    net_file,
    net_df=None,
    anno=None,
    MI_thres=0,
    regul_size=50,
    normalize_MI_per_regulon=True
):
    """\
    Process an output from ARACNe3 to return a pd.DataFrame describing a gene
    regulatory network with suitable columns for conversion to an object of the
    Interactome class.

    Parameters
    ----------
    net_file
        A string containing the path to the ARACNe3 output
    net_df (default: None)
        Whether to passt a pd.DataFrame instead of the path
    anno (default: None)
        Gene ID annotation
    MI_thres (default: 0)
        Threshold on Mutual Information (MI) to select the regulators and target pairs
    regul_size (default: 50)
        Number of (top) targets to include in each regulon
    normalize_MI_per_regulon (default: True)
        Whether to normalize MI values each regulon by the maximum value

    Returns
    -------
    A pd.DataFrame containing an ARACNe3-inferred gene regulatory network with the
    following 4 columns: "regulator", "target", "mor" (mode of regulation) and "likelihood".
    """
    return _aracne3_to_regulon(
        net_file,
        net_df,
        anno,
        MI_thres,
        regul_size,
        normalize_MI_per_regulon
    )

# def nes_to_neg_log(adata, layer = None, key_added = None):
#     """\
#     Transform VIPER-computed NES into -log10(p-value).
#
#     Parameters
#     ----------
#     adata
#         Gene expression, protein activity or pathways stored in an anndata
#         object, or a pandas dataframe containing input data, where rows are
#         observations/samples (e.g. cells or groups) and columns are features
#         (e.g. proteins or pathways).
#     layer : (default: None)
#         Entry of layers to tranform.
#     key_added : (default: None)
#         Name of layer to save result in a new layer instead of adata.X.
#
#     Returns
#     -------
#     Saves the input data as a transformed version. If key_added is specified,
#     saves the results in adata.layers[key_added].
#     """
#     _nes_to_neg_log(adata, layer, key_added)

def nes_to_pval(
    adata,
    layer = None,
    key_added = None,
    lower_tail=True,
    adjust=True,
    axs=1,
    neg_log = False,
    copy = False
):
    """\
    Transform VIPER-computed NES into p-values.

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata
        object, or a pandas dataframe containing input data, where rows are
        observations/samples (e.g. cells or groups) and columns are features
        (e.g. proteins or pathways).
    layer : (default: None)
        Entry of layers to tranform.
    key_added : (default: None)
        Name of layer to save result in a new layer instead of adata.X.
    lower_tail: default (True)
    	If `True` (default), probabilities are P(X <= x)
    	If `False`, probabilities are P(X > x)
    adjust (default: True)
        If `True`, returns adjusted p values using FDR Benjamini-Hochberg procedure.
        If `False`, does not adjust p values
    axs (default: 1)
        axis along which to perform the p-value correction (Used only if the input is a pd.DataFrame).
        Possible values are 0 or 1.
    neg_log : (default: False)
        Whether to transform VIPER-computed NES into -log10(p-value).
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    Saves the input data as a transformed version. If key_added is specified,
    saves the results in adata.layers[key_added].
    """
    return _nes_to_pval(adata, layer, key_added, lower_tail, adjust, axs, neg_log, copy)

def repr_subsample(adata,
                   pca_slot="X_pca",
                   size=1000,
                   seed=0,
                   key_added = "repr_subsample",
                   eliminate = False,
                   verbose=True,
                   njobs=1,
                   copy = False):
    """\
    A tool for create a subsample of the input data such it is well
    representative of all the populations within the input data rather than
    being a random sample. This is accomplished by pairing samples together in
    an iterative fashion until the desired sample size is reached.

    Parameters
    ----------
    adata
        An anndata object containing a distance object in adata.obsp.
    pca_slot (default: "X_pca")
        The slot in adata.obsm where the PCA object is stored. One way of
        generating this object is with sc.pp.pca.
    size (default: 1000)
        The size of the representative subsample
    eliminate (default: False)
        Whether to trim down adata to the subsample (True) or leave the
        subsample as an annotation in adata.obs[key_added].
    seed (default: 0)
        The random seed used when taking samples of the data.
    verbose (default: True)
        Whether to provide runtime information.
    njobs (default: 1)
        The number of cores to use for the analysis. Using more than 1 core
        (multicore) speeds up the analysis.
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    When copy is False, saves the subsample annotation in adata.var[key_added].
    When copy is True, return an anndata with this annotation. When eliminate is
    True, modify the adata by subsetting it down to the subsample.
    """
    return _representative_subsample_anndata(
        adata,
        pca_slot,
        size,
        exact_size = True,
        seed = seed,
        key_added = key_added,
        eliminate = eliminate,
        verbose = verbose,
        njobs = njobs,
        copy = copy
    )


# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ----------------------------- ** METACELL FUNC ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def repr_metacells(
    adata,
    counts = None,
    pca_slot = "X_pca",
    dist_slot = "corr_dist",
    clusters_slot = None,
    score_slot = None,
    score_min_thresh = None,
    size = 500,
    n_cells_per_metacell = None,
    min_median_depth = 10000,
    perc_data_to_use = None,
    perc_incl_data_reused = None,
    seed = 0,
    key_added = "metacells",
    verbose = True,
    njobs = 1,
    copy = False
):
    """\
    A tool for create a representative selection of metacells from the data that
    aims to maximize reusing samples from the data, while simultaneously
    ensuring that all neighbors are close to the metacell they construct.
    When using this function, exactly two of the following parameters must be
    set: size, min_median_depth or n_cells_per_metacell, perc_data_to_use or
    perc_incl_data_reused.
    Note that min_median_depth and n_cells_per_metacell cannot both be set
    at the same time, since they directly relate (e.g. higher n_cells_per_metacell
    means more neighbors are used to construct a single metacell, meaning each
    metacell will have more counts, resulting in a higher median depth).
    Note that perc_data_to_use and perc_incl_data_reused cannot both be set
    at the same time, since they directly relate (e.g. higher perc_data_to_use
    means you include more data, which means it's more likely to reuse more
    data, resulting in a higher perc_incl_data_reused).

    Parameters
    ----------
    adata
        An anndata object containing a distance object in adata.obsp.
    counts (default: None)
        A pandas DataFrame or AnnData object of unnormalized gene expression
        counts that has the same samples in the same order as that of adata.
        If counts are left as None, adata must have counts stored in adata.raw.
    pca_slot (default: "X_pca")
        The slot in adata.obsm where the PCA object is stored. One way of
        generating this object is with sc.pp.pca.
    dist_slot (default: "corr_dist")
        The slot in adata.obsp where the distance object is stored. One way of
        generating this object is with pyviper.pp.corr_distance.
    clusters_slot (default: None)
        The slot in adata.obs where cluster labels are stored. Cluster-specific
        metacells will be generated using the same parameters with the results
        for each cluster being stored separately in adata.uns.
    score_slot (default: None)
        The slot in adata.obs where a score used to determine and filter cell
        quality are stored (e.g. silhouette score).
    score_min_thresh (default: None)
        The score from adata.obs[score_slot] that a cell must have at minimum to
        be used for metacell construction (e.g. 0.25 is the rule of thumb for
        silhouette score).
    size (default: None)
        A specific number of metacells to generate. If left as None,
        perc_data_to_use or perc_incl_data_reused can be used to specify the size
        when n_cells_per_metacell or min_median_depth is given.
    n_cells_per_metacell (default: None)
        The number of cells that should be used to generate single metacell.
        Note that this parameter and min_median_depth cannot both be set as
        they directly relate: e.g. higher n_cells_per_metacell leads to higher
        min_median_depth. If left as None, perc_data_to_use or
        perc_incl_data_reused can be used to specify n_cells_per_metacell when
        size is given.
    min_median_depth (default: 10000)
        The desired minimum median depth for the metacells (indirectly specifies
        n_cells_per_metacell). The default is set to 10000 as this is recommend
        by PISCES[1]. Note that this parameter and n_cells_per_metacell cannot
        both be set as they directly relate: e.g. higher min_median_depth leads
        to higher n_cells_per_metacell.
    perc_data_to_use (default: None)
        The percent of the total amount of provided samples that will be used in
        the creation of metacells. Note that this parameter and
        perc_incl_data_reused cannot both be set as they directly relate: e.g.
        higher perc_data_to_use leads to higher perc_incl_data_reused.
    perc_incl_data_reused (default: None)
        The percent of samples that are included in the creation of metacells
        that will be reused (i.e. used in more than one metacell). Note that this
        parameter and perc_data_to_use cannot both be set as they directly relate:
        e.g. higher perc_incl_data_reused leads to higher perc_data_to_use.
    seed (default: 0)
        The random seed used when taking samples of the data.
    key_added (default: "metacells")
        The name of the slot in the adata.uns to store the output.
    verbose (default: True)
        Whether to provide runtime information and quality statistics.
    njobs (default: 1)
        The number of cores to use for the analysis. Using more than 1 core
        (multicore) speeds up the analysis.
    copy (default: False)
        Determines whether a copy of the input AnnData is returned.

    Returns
    -------
    Saves the metacells as a pandas dataframe in adata.uns[key_added]. Attributes
    that contain parameters for and statistics about the construction of the
    metacells are stored in adata.uns[key_added].attrs. Set copy = True to
    return a new AnnData object.

    Citations
    -------
    Obradovic, A., Vlahos, L., Laise, P., Worley, J., Tan, X., Wang, A., &
    Califano, A. (2021). PISCES: A pipeline for the systematic, protein activity
    -based analysis of single cell RNA sequencing data. bioRxiv, 6, 22.
    """

    # Here is a summary of parameter choices and their affects:
	# size + min_median_depth
	# 	--> calculate n_cells_per_metacell
	# size + n_cells_per_metacell
	# 	--> N/A
	# perc_data_used + size
	# 	--> optimize n_cells_per_metacell
	# perc_data_used + n_cells_per_metacell
	# 	--> optimize size
	# perc_data_used + min_median_depth
	# 	--> calculate n_cells_per_metacell
	# 	--> optimize size
	# max_perc_included_reused + size
	# 	--> optimize n_cells_per_metacell
	# max_perc_included_reused + n_cells_per_metacell
	# 	--> optimize size
	# max_perc_included_reused + min_median_depth
	# 	--> calculate n_cells_per_metacell
	# 	--> optimize size

    # Removed single parameter choices:
    # size
	# 	--> assume n_cells_per_metacell = total_samples/size
	# n_cells_per_metacell
	# 	--> assume size = total_samples/n_cells_per_metacell
	# min_median_depth
	# 	--> calculate n_cells_per_metacell
	# 	--> assume size = total_samples/n_cells_per_metacell

    if score_slot in adata.obs.columns is not None:
        adata = adata[adata.obs[score_slot].values >= score_min_thresh,:].copy()

    return _representative_metacells_multiclusters(
        adata,
        counts,
        pca_slot,
        dist_slot,
        clusters_slot,
        size,
        n_cells_per_metacell,
        min_median_depth,
        perc_data_to_use,
        perc_incl_data_reused,
        exact_size = True,
        seed = seed,
        key_added = key_added,
        verbose = verbose,
        njobs = njobs,
        smart_sample = True,
        copy = copy
    )
