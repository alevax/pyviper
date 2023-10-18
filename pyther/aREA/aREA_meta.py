### ---------- IMPORT DEPENDENCIES ----------
import numpy as np
import pandas as pd
from .._helpers_meta import get_resized_mats
from .aREA_classic import aREA_classic
from ..interactome import Interactome

### ---------- EXPORT LIST ----------
__all__ = ['aREA']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def consolidate_meta_aREA_results(netMets, mvws = 1, verbose = True):
    if(verbose): print("Integrating NES matrices together with mvws=" + str(mvws) + "...")

    # Resize the matrices so they all share the same shape, row names and column names
    resized_nes_list = get_resized_mats(netMets, empty_value = 0)

    # Stack the data now that it's the same size
    nes_arrays_list = [df.values for df in resized_nes_list]
    stacked_nes = np.stack(nes_arrays_list, axis=-1)

    # Compute the weights - creating a new set of stacks
    ws = np.abs(stacked_nes)**mvws

    # Multiply the stacks together; sum across the stacks; divide the sums
    # Note that with mvws=1, multiplying weights and dividing here cancel out
    meta_array = np.sum(stacked_nes*ws, axis=-1)/np.sum(ws, axis=-1)

    # 0 values that appears a NES scores will create some weights that are
    #  also 0. Thus, all NaN values that appear are 0/0 division that should
    # be corrected to 0
    meta_array = np.nan_to_num(meta_array)

    # Put everything in a nice user-friendly DataFrame
    preOp = pd.DataFrame(meta_array,
                         index=resized_nes_list[0].index,
                         columns=resized_nes_list[0].columns).sort_index()
    return preOp

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def aREA(gex_data, interactome, layer = None, eset_filter = False, min_targets = 30, mvws = 1, verbose = True):
    """\
    Allows the individual to infer normalized enrichment scores from gene
    expression data using the Analytical Ranked Enrichment Analysis (aREA)[1]
    function.

    It is the original basis of the VIPER (Virtual Inference of Protein-activity
    by Enriched Regulon analysis) algorithm.

    The Interactome object must not contain any targets that are not in the
    features of gex_data. This can be accomplished by running:
        interactome.filter_targets(gex_data.var_names)
    It is highly recommended to do this on the unPruned network and then prune to
    ensure the pruned network contains a consistent number of targets per
    regulator, all of which exist within gex_data. A consistent number of
    targets allows regulators to have NES scores that are comparable to one
    another. A regulator that has more targets than others will have "boosted"
    NES scores, such that they cannot be compared to those with fewer targets.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy).
    interactome
        An object of class Interactome.
    layer (default: None)
        The layer in the anndata object to use as the gene expression input.
    eset_filter (default: False)
        Whether to filter out genes not present in the interactome (True) or to
        keep this biological context (False). This will affect gene rankings.
    min_targets (default: 30)
        The minimum number of targets that each regulator in the interactome
        should contain. Regulators that contain fewer targets than this minimum
        will be culled from the network (via the Interactome.cull method). The
        reason users may choose to use this threshold is because adequate
        targets are needed to accurately predict enrichment.
    mvws (default: 1)
        Number indicating either the exponent score for the metaViper weights.
        Roughly, a lower number (e.g. 1) results in networks being treated as
        a consensus network (useful for multiple networks of the same celltype
        with the same epigenetics), while a higher number (e.g. 10) results in
        networks being treated as separate (useful for multiple networks of
        different celltypes with different epigenetics).
    verbose (default: True)
        Whether extended output about the progress of the algorithm should be
        given.

    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.

    Citations
    -------
    [1] Alvarez, M. J., Shen, Y., Giorgi, F. M., Lachmann, A., Ding, B. B., Ye,
    B. H., & Califano, A. (2016). Functional characterization of somatic
    mutations in cancer using network-based inference of protein activity.
    Nature genetics, 48(8), 838-847.
    """

    if not isinstance(mvws, int): raise ValueError("mvws is not of type int.")
    # We want all if/else conditions in case
    # users or testers run this function directly

    if isinstance(interactome, Interactome):
        preOp = aREA_classic(gex_data, interactome, layer, eset_filter, min_targets, verbose)
    elif len(interactome) == 1:
        preOp = aREA_classic(gex_data, interactome[0], layer, eset_filter, min_targets, verbose)
    else:
        # netMets = [aREA(gex_data, iObj, eset_filter, layer) for iObj in interactome]
        netMets = []
        tot_nets = len(interactome)
        n_completed_nets = 0
        if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
        for iObj in interactome:
            netMets.append(aREA_classic(gex_data, iObj, layer, eset_filter, min_targets, verbose))
            n_completed_nets = n_completed_nets + 1
            if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
        preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)

    return preOp
