### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
from .NaRnEA_classic import *
from .._helpers_meta import *
from ..interactome import Interactome


### ---------- EXPORT LIST ----------
__all__ = ['NaRnEA']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def get_pes_list(results):
    # Iterate through each DataFrame to collect unique gene names
    pes_list = []
    for res in results:
        pes_list.append(res['pes'])
    return pes_list

def get_all_regs_in_NaRnEA_list(results):
    pes_list = get_pes_list(results)
    return get_all_regs(pes_list)

def get_resized_pes(results):
    # Get list of PES dataframes
    pes_list = get_pes_list(results)
    resized_pes_list = get_resized_mats(pes_list, empty_value = np.nan)
    return resized_pes_list

def find_max_absolute_value_df(dataframes):
    # Stack the DataFrames into a 3D NumPy array
    stacked_arrays = [df.values for df in dataframes]
    stacked_data = np.stack(stacked_arrays, axis=-1)

    # Calculate the absolute maximum along the last axis (axis=-1)
    absolute_max_indices = np.nanargmax(np.abs(stacked_data), axis=-1).astype(float)

    # Calculate the existance of any nan values along the last axis (axis=-1)
    nan_present_indices = np.isnan(stacked_data).any(axis=-1)

    # Replace values in absolute_max_indices with NaN wherever NaN is present in
    # nan_present_indices to remove voting power for those regulators
    absolute_max_indices[nan_present_indices] = np.nan

    # Create a DataFrame with the same row and column names
    result_df = pd.DataFrame(absolute_max_indices, index=dataframes[0].index, columns=dataframes[0].columns)

    return result_df

def calculate_value_proportions(result_df, n_nets):
    # Convert the DataFrame to a NumPy array
    result_array = result_df.values

    # Initialize an empty NumPy array to store proportions
    proportions_array = np.zeros((result_array.shape[0], n_nets))

    # Count the occurrences of each value in each row, ignore NaN values
    for i, row in enumerate(result_array):
        integer_values = row[np.isfinite(row) & (row == np.floor(row))]
        unique_values, counts = np.unique(integer_values, return_counts=True)
        proportions_array[i, unique_values.astype(int)] = counts / len(integer_values)

    # Convert the proportions array to a DataFrame
    proportions_df = pd.DataFrame(proportions_array, columns=range(n_nets), index=result_df.index)

    return proportions_df

def get_net_weight(results):
    # Get the total number of networks
    n_nets = len(results)
    # Resize the PES matrices so they have the same size, column names and row names
    resized_pes_list = get_resized_pes(results)
    # Stack the resized PES matrices: along the z axis, calculate position of max abs value
    max_abs_vals_df = find_max_absolute_value_df(resized_pes_list)
    # For each sample, calculate the proportion of max abs values from each network
    net_weight = calculate_value_proportions(max_abs_vals_df, n_nets).sort_index()

    net_weight.index.name = 'index'
    net_weight.columns.name = 'net'

    return net_weight

def integrate_NaRnEA_xes_mats(results, bg_matrix, net_weight, xes_type = 'pes'):

    pre_xes = bg_matrix.copy()
    xes_dom = bg_matrix.copy()
    n_regs = bg_matrix.shape[1]

    for i in range(0,len(results)):
        all_regs_xes = bg_matrix.copy() + results[i][xes_type]
        all_regs_xes.fillna(0, inplace=True)
        all_regs_xes.sort_index(inplace=True)
        if xes_type == 'nes':
            xes_dom = xes_dom + ( (all_regs_xes!=0) * net_weight[i].to_numpy()[:, np.newaxis] * np.ones((1, n_regs)) )**2
        else: #xes_type == 'pes':
            xes_dom = xes_dom + ( (all_regs_xes!=0) * net_weight[i].to_numpy()[:, np.newaxis] * np.ones((1, n_regs)) )
        net_weight_array = net_weight[i].to_numpy()[:, np.newaxis] * np.ones((1, n_regs))
        pre_xes =  pre_xes + all_regs_xes * net_weight_array

    if xes_type == 'nes':
        pre_xes = pre_xes/np.sqrt(xes_dom)
    else: #xes_type == 'pes':
        pre_xes = pre_xes/xes_dom

    # There must be regulator activity in some cells where neither network provides nonzero value.
    # Imagine we have a regulator A. A has a regulon A_1 and A_2 in net_1 and net_2 respectively.
    # Imagine, the targets of A_1 and the targets of A_2 are both 0 in samples S1, S2, ...., Sn.
    # For these samples pre_xes will have 0 and these 0s will be present in xes_dom in the same place.
    # This is because we start with bg_matrix (0 matrix) and add to it wherever the xes matrices are not 0 (all_regs_xes!= 0)
    # so if the xes matrices have 0s in them, this will result in 0s remaining in the dom_matrix.
    # Hence fillna fixes 0/0, which should just remain at 0.
    # The following print statements will show the same numbers of 0s.
    # print(np.count_nonzero(pre_xes == 0))
    # print(np.count_nonzero(xes_dom == 0))
    xes = pre_xes.fillna(0)

    return(xes)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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
    verbose (default: True)
        Whether extended output about the progress of the algorithm should be
        given.

    Returns
    -------
    A dictionary containing :class:`~numpy.ndarray` containing NES values
    (key: 'nes') and PES values (key: 'pes').

    Citations
    -------
    [1] Griffin, A. T., Vlahos, L. J., Chiuzan, C., & Califano, A. (2023). NaRnEA:
    An Information Theoretic Framework for Gene Set Analysis. Entropy, 25(3), 542.
    """
    pd.options.mode.chained_assignment = None

    if isinstance(interactome, Interactome):
        return NaRnEA_classic(gex_data, interactome, layer, eset_filter, min_targets, verbose)
    elif len(interactome) == 1:
        return NaRnEA_classic(gex_data, interactome[0], layer, eset_filter, min_targets, verbose)
    else:
        # results = [NaRnEA(gex_data, iObj, verbose = verbose) for iObj in interactome]
        results = []
        tot_nets = len(interactome)
        n_completed_nets = 0
        if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
        for iObj in interactome:
            results.append(NaRnEA_classic(gex_data, iObj, layer, eset_filter, min_targets, verbose))
            n_completed_nets = n_completed_nets + 1
            if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")

    if verbose: print('Integrating results')

    net_weight = get_net_weight(results)
    all_regs = get_all_regs_in_NaRnEA_list(results)
    bg_matrix = pd.DataFrame(0,index = results[0]['nes'].index, columns = all_regs)

    nes = integrate_NaRnEA_xes_mats(results, bg_matrix, net_weight, xes_type = 'nes')
    pes = integrate_NaRnEA_xes_mats(results, bg_matrix, net_weight, xes_type = 'pes')

    result = {"nes": nes, "pes": pes}
    return result
