### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
from ._helpers import *
from .aREA.aREA_meta import aREA
from .NaRnEA.NaRnEA_meta import NaRnEA

### ---------- EXPORT LIST ----------
__all__ = ['vithon']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def sample_ttest(i,array):
    return ttest_1samp((array[i] - np.delete(array, i, 0)), 0).statistic

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def vithon(gex_data,
           interactome,
           layer=None,
           eset_filter=True,
           method=None,  # [None, "scale", "rank", "mad", "ttest"],
           enrichment='area',  # [None, 'area','narnea'],
           mvws=1,
           min_targets=30,
           njobs=1,
           verbose=True,
           output_as_anndata=True,
           transfer_obs=True,
           store_input_data=True
           ):
    """\
    Allows the individual to infer normalized enrichment scores and proportional
    enrichment scores from gene expression data using the Nonparametric
    analytical rank-based enrichment analysis (NaRnEA) function.

    NaRnEA is an updated basis for the VIPER (Virtual Inference of Protein-activity
    by Enriched Regulon analysis) algorithm.

    Griffin, A. T., Vlahos, L. J., Chiuzan, C., & Califano, A. (2023). NaRnEA:
    An Information Theoretic Framework for Gene Set Analysis. Entropy, 25(3), 542.

    The Interactome object must not contain any targets that are not in the
    features of gex_data. This can be accomplished by running:
        interactome.filter_targets(gex_data.var_names)
    It is highly recommend to do this on the unPruned network and then prune to
    ensure the pruned network contains a consistent number of targets per
    regulator, allow of which exist within gex_data.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy).
    interactome
        An object of class Interactome or a list of Interactome objects.
    layer (default: None)
        The layer in the anndata object to use as the gene expression input.
    eset_filter (default: False)
        Whether to filter out genes not present in the interactome (True) or to
        keep this biological context (False). This will affect gene rankings.
    method (default: None)
        A method used to create a gene expression signature from gex_data.X. The
        default of None is used when gex_data.X is already a gene expression
        signature. Alternative inputs include "scale", "rank", "mad", and "ttest".
    enrichment (default: 'area')
        The algorithm to use to calculate the enrichment. Choose betweeen
        Analytical Ranked Enrichment Analysis (aREA) and Nonparametric
        Analytical Rank-based Enrichment Analysis (NaRnEA) function. Default ='area',
        alternative = 'narnea'.
    mvws (default: 1)
        Number indicating either the exponent score for the metaViper weights.
        These are only applicable when enrichment = 'area' and are not used when
        enrichment = 'narnea'. Roughly, a lower number (e.g. 1) results in
        networks being treated as a consensus network (useful for multiple
        networks of the same celltype with the same epigenetics), while a higher
        number (e.g. 10) results in networks being treated as separate (useful
        for multiple networks of different celltypes with different epigenetics).
    min_targets (default: 30)
        The minimum number of targets that each regulator in the interactome
        should contain. Regulators that contain fewer targets than this minimum
        will be culled from the network (via the Interactome.cull method). The
        reason users may choose to use this threshold is because adequate
        targets are needed to accurately predict enrichment.
    njobs (default: 1)
        For metaViper. When using multiple networks, assign jobs to different
        networks.
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
    A dictionary containing :class:`~numpy.ndarray` containing NES values
    (key: 'nes') and PES values (key: 'pes') when output_as_anndata==False and
    enrichment = "narnea".
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values
    when output_as_anndata==False and enrichment = "area".
    An anndata object containin NES values in .X. Will contain PES values in the
    layer 'pes' when enrichment = 'narnea'. Will contain .gex_data and/or
    .pax_data in the unstructured data slot (.uns) when store_input_data = True.
    Will contain identical .obs to the input anndata when transfer_obs = True.
    """

    gex_data_original = gex_data
    gex_data = gex_data_original.copy()

    pd.options.mode.chained_assignment = None

    if verbose: print("Preparing the association scores")

    ''' not in current version 3.8 of python , only in python 3.10
    match method:

        case 'scale': #scale each gene in all sample
            gex_data.X = (gex_data.X - np.mean(gex_data.X,axis=0))/np.std(gex_data.X,axis=0)
        case 'rank':
            gex_data.X = rankdata(gex_data.X,axis=0)*(np.random.random(gex_data.X.shape)*2/10-0.1)
        case 'mad':
            median = np.median(gex_data.X,axis=0)
            gex_data.X = (gex_data.X-median)/np.median(np.abs(gex_data.X-median))
        case 'ttest':
            # I dont quite understand this part maybe ask later:
            gex_data.X = gex_data.X
    '''
    if method == 'scale':
        gex_data.X = (gex_data.X - np.mean(gex_data.X,axis=0))/np.std(gex_data.X,axis=0)
    elif method == 'rank':
        gex_data.X = rankdata(gex_data.X,axis=0)*(np.random.random(gex_data.X.shape)*2/10-0.1)
    elif method == 'mad':
        median = np.median(gex_data.X,axis=0)
        gex_data.X = (gex_data.X-median)/np.median(np.abs(gex_data.X-median))
    elif method == 'ttest':
        gex_data.X = np.array([sample_ttest(i, gex_data.X.copy()) for i in range(gex_data.shape[0])])

    if enrichment is None:
        enrichment = 'narnea'
    else:
        enrichment = enrichment.lower()
    if enrichment == 'area':
        if verbose: print("Computing regulons enrichment with aREA")
        preOp = aREA(gex_data, interactome, layer, eset_filter,
                     min_targets, mvws, njobs, verbose)
    elif enrichment == 'narnea':
        if verbose: print("Computing regulons enrichment with NaRnEa")
        preOp = NaRnEA(gex_data, interactome, layer, eset_filter,
                       min_targets, njobs, verbose)
    else:
        raise ValueError("Unsupported enrichment type:" + str(enrichment))

    if output_as_anndata == False:
        op = preOp
    else: #output_as_anndata == True:
        if enrichment == 'area':
            op = mat_to_anndata(preOp)
        else: #enrichment == 'narnea':
            op = mat_to_anndata(preOp["nes"])
            op.layers['pes'] = preOp["pes"]

        if transfer_obs is True:
            op.obs = op.obs.join(gex_data.obs)
        if store_input_data is True:
            # If input data was pax_data for pathway enrichment
            if 'gex_data' in gex_data_original.uns:
                op.uns['gex_data'] = gex_data_original.uns['gex_data']
                op.uns['pax_data'] = gex_data_original
            else:
                op.uns['gex_data'] = gex_data_original
    return op #final result
