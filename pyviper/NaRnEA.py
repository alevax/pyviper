from ._NaRnEA.NaRnEA_meta import NaRnEA_meta

def NaRnEA(gex_data, interactome, layer = None, eset_filter = False, min_targets = 30, verbose = True):
    """\
    Allows the individual to infer normalized enrichment scores and proportional
    enrichment scores from gene expression data using the Nonparametric
    Analytical Rank-based Enrichment Analysis (NaRnEA)[1] function.

    NaRnEA is an updated basis for the VIPER (Virtual Inference of
    Protein-activity by Enriched Regulon analysis) algorithm.

    The Interactome object must not contain any targets that are not in the
    features of gex_data. This can be accomplished by running:
        interactome.filter_targets(gex_data.var_names)
    It is highly recommend to do this on the unPruned network and then prune to
    ensure the pruned network contains a consistent number of targets per
    regulator, all of which exist within gex_data. A regulator that has more
    targets than others will have "boosted" NES scores, such that they cannot be
    compared to those with fewer targets.

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
    verbose : default: True
        Whether extended output about the progress of the algorithm should be
        given.

    Returns
    -------
    A dictionary containing :class:`~numpy.ndarray` containing NES values (key: 'nes') and PES values (key: 'pes').

    References
    ----------
    [1] Griffin, A. T., Vlahos, L. J., Chiuzan, C., & Califano, A. (2023). NaRnEA: An Information Theoretic Framework for Gene Set Analysis. Entropy, 25(3), 542.
    """
    return NaRnEA_meta(
        gex_data,
        interactome,
        layer,
        eset_filter,
        min_targets,
        verbose
    )
