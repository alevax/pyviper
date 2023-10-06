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

def aREA(gex_data, interactome, eset_filter = False, layer = None, min_targets = 30, mvws = 1, njobs = 1, verbose = True):
    if not isinstance(mvws, int): raise ValueError("mvws is not of type int.")
    # We want all if/else conditions in case
    # users or testers run this function directly

    if isinstance(interactome, Interactome):
        preOp = aREA_classic(gex_data, interactome, eset_filter, layer, min_targets, verbose)
    elif len(interactome) == 1:
       preOp = aREA_classic(gex_data, interactome[0], eset_filter, layer, min_targets, verbose)
    elif njobs == 1:
        # netMets = [aREA(gex_data, iObj, eset_filter, layer) for iObj in interactome]
        netMets = []
        tot_nets = len(interactome)
        n_completed_nets = 0
        if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
        for iObj in interactome:
          netMets.append(aREA_classic(gex_data, iObj, eset_filter, layer, min_targets, verbose))
          n_completed_nets = n_completed_nets + 1
          if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
        preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)
    else:
        joblib_verbose = 0
        if verbose:
            print("Computing regulons enrichment with aREA")
            joblib_verbose = 11
        # n_jobs need to be decided.
        netMets = Parallel(n_jobs = njobs, verbose = joblib_verbose)(
            (delayed)(aREA)(gex_data, iObj, eset_filter, layer, verbose)
            for iObj in interactome
            )
        preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)
    return preOp
