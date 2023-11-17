### ---------- IMPORT DEPENDENCIES ----------
import numpy as np
import pandas as pd
from .._helpers_meta import get_resized_mats
from .aREA_classic import aREA_classic
from ..interactome import Interactome
from scipy.stats import rankdata

### ---------- EXPORT LIST ----------
__all__ = ['aREA']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def consolidate_meta_aREA_results_by_weights(netMats, mvws = 1):
    # Resize the matrices so they all share the same shape, row names and column names
    resized_nes_list = get_resized_mats(netMats, empty_value = 0)

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

# This code is for get_net_index_assignments when mvws = "auto"
def __rank_column(column): return rankdata(column, method='ordinal')
def __get_net_score_per_sample(vpdf, n_mrs = 50):
    vpmat = vpdf.values
    n_cols = vpmat.shape[1]

    # For each sample, sort MRs by distance.
    vpmat_ranked = np.apply_along_axis(
        __rank_column,
        axis=1,
        arr=np.array(vpmat)
    )

    # Sort each row of the input array and get the indices of the sorted elements
    sorted_indices = np.argsort(vpmat_ranked, axis=1)

    # Slice the sorted indices to keep only the top knn indices for each row
    top_n_indices = sorted_indices[:, (n_cols-n_mrs):n_cols]

    # Get the NES values for each of those top_n_indices
    top_n_NES = np.take_along_axis(vpmat, top_n_indices, axis=1)

    # Take the sum of those top_n_indices
    top_n_NES_means = np.mean(top_n_NES, axis=1)

    # Do the same but for the bottom candidate MRs
    bottom_n_indices = sorted_indices[:, :n_mrs]
    bottom_n_NES = np.take_along_axis(vpmat, bottom_n_indices, axis=1)
    bottom_n_NES_means = np.mean(bottom_n_NES, axis=1)

    net_scores = (top_n_NES_means + bottom_n_NES_means*-1)/2
    return net_scores
def __get_net_index_assignments(netMats):
    n_nets = len(netMats)
    nets_scores_array = np.zeros((netMats[0].shape[0], n_nets))
    for i in range(n_nets):
        nets_scores_array[:,i] = __get_net_score_per_sample(netMats[i], n_mrs = 50)
    net_index_assignments = max_indices = np.argmax(nets_scores_array, axis=1)
    return net_index_assignments
def consolidate_meta_aREA_results_by_network_matching(netMats):
    # Assign each sample to a network by the top candidate MR scores
    net_index_assignments = __get_net_index_assignments(netMats)

    # Resize the matrices so they all share the same shape, row names and column names
    resized_nes_list = get_resized_mats(netMats, empty_value = 0)

    # Stack the arrays
    stacked_arrays = np.stack(resized_nes_list)

    # Index the stack by the net_index_assignments to select rows based on matching samples ot networks
    meta_array = stacked_arrays[net_index_assignments, np.arange(len(net_index_assignments)), :]

    # Save the final array in the pd.DataFrame format
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

    # We want all if/else conditions in case
    # users or testers run this function directly
    if isinstance(mvws, str) and mvws != "auto":
        # network manual assignment
        if(verbose): print("Running VIPER using manual network matching annotation...")
        netMats = []
        tot_nets = len(interactome)
        n_completed_nets = 0
        net_assignments = gex_data.obs[mvws].values
        if isinstance(net_assignments[0], str):
            # Create net_index_assignments by replacing each on with the appropriate network names
            int_names_in_index_order = [iObj.name for iObj in interactome]
            name_positions = {name: position for position, name in enumerate(int_names_in_index_order)}
            # Use vectorized operations to get positions for each item in net_assignments
            net_index_assignments = np.vectorize(name_positions.get)(net_assignments)
        elif isinstance(net_assignments[0], np.integer):
            # Treat these as index_assignments
            net_index_assignments = net_assignments
        else:
            raise ValueError("Manual network assignments must correspond to network names\n" +
                             "(str) or be indices to the input list interactome (int)")
        for i in range(tot_nets):
            iObj = interactome[i]
            gex_data_subset = gex_data[net_index_assignments==i,]
            if gex_data_subset.shape[0] == 0: continue # for multi-core
            netMats.append(aREA_classic(gex_data_subset, iObj, layer, eset_filter, min_targets, verbose))
            n_completed_nets = n_completed_nets + 1
            if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")

        resized_nes_list = get_resized_mats(netMats, empty_value = np.nan)
        preOp = pd.concat(resized_nes_list, ignore_index=True)
    else:
        if not (isinstance(mvws, int) or mvws == "auto"):
            raise ValueError('mvws must be "auto", an int, or a vector of manual assignments.')

        if isinstance(interactome, Interactome):
            preOp = aREA_classic(gex_data, interactome, layer, eset_filter, min_targets, verbose)
        elif len(interactome) == 1:
            preOp = aREA_classic(gex_data, interactome[0], layer, eset_filter, min_targets, verbose)
        else:
            # netMats = [aREA(gex_data, iObj, eset_filter, layer) for iObj in interactome]
            netMats = []
            tot_nets = len(interactome)
            n_completed_nets = 0
            if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
            for iObj in interactome:
                netMats.append(aREA_classic(gex_data, iObj, layer, eset_filter, min_targets, verbose))
                n_completed_nets = n_completed_nets + 1
                if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
            # Consolidate the data together
            if isinstance(mvws, int):
                if(verbose): print("Integrating NES matrices together with mvws=" + str(mvws) + "...")
                preOp = consolidate_meta_aREA_results_by_weights(netMats, mvws)
            elif mvws == "auto":
                if(verbose): print("Integrating NES matrices together by automatic network matching...")
                preOp = consolidate_meta_aREA_results_by_network_matching(netMats)

    return preOp
