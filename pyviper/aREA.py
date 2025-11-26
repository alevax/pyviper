from ._aREA.aREA_meta import aREA_meta

def aREA(gex_data, interactome, layer = None, eset_filter = False, min_targets = 30, mvws = 1, verbose = True):
    """
    Allows the individual to infer normalized enrichment scores from gene
    expression data using the Analytical Ranked Enrichment Analysis (aREA)[1]
    function.
    
    It is the original basis of the VIPER (Virtual Inference of Protein-activity
    by Enriched Regulon analysis) algorithm.

    The Interactome object must not contain any targets that are not in the
    features of gex_data. This can be accomplished by running:

        interactome.filter_targets(gex_data.var_names)
    
    It is highly recommended to do this on the unPruned network and then prune
    to ensure the pruned network contains a consistent number of targets per
    regulator, all of which exist within gex_data. A consistent number of
    targets allows regulators to have NES scores that are comparable to one
    another. A regulator that has more targets than others will have "boosted"
    NES scores, such that they cannot be compared to those with fewer targets.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy) or in a
        pd.DataFrame.
    interactome
        An object of class Interactome.
    layer : default: None
        The layer in the anndata object to use as the gene expression input.
    eset_filter : default: False
        Whether to filter out genes not present in the interactome (True) or to
        keep this biological context (False). This will affect gene rankings.
    min_targets : default: 30
        The minimum number of targets that each regulator in the interactome
        should contain. Regulators that contain fewer targets than this minimum
        will be culled from the network (via the Interactome.cull method). The
        reason users may choose to use this threshold is because adequate
        targets are needed to accurately predict enrichment.
    mvws : default: 1
        (A) Number indicating either the exponent score for the metaViper weights.
        These are only applicable when enrichment = 'area' and are not used when
        enrichment = 'narnea'. Roughly, a lower number (e.g. 1) results in
        networks being treated as a consensus network (useful for multiple
        networks of the same celltype with the same epigenetics), while a higher
        number (e.g. 10) results in networks being treated as separate (useful
        for multiple networks of different celltypes with different epigenetics).
        (B) The name of a column in gex_data that contains the manual assignments
        of samples to networks using list position or network names.
        (C) "auto": assign samples to networks based on how well each
        network allows for sample enrichment.
    verbose : default: True
        Whether extended output about the progress of the algorithm should be
        given.

    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.

    References
    ----------
    [1] Alvarez, M. J., Shen, Y., Giorgi, F. M., Lachmann, A., Ding, B. B., Ye, B. H., & Califano, A. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. Nature genetics, 48(8), 838-847.
    """
    return aREA_meta(
        gex_data,
        interactome,
        layer,
        eset_filter,
        min_targets,
        mvws,
        verbose
    )
