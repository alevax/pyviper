### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from scipy.special import ndtri # from scipy.stats import norm
import warnings
import time
import torch
from anndata import AnnData


### ---------- EXPORT LIST ----------
__all__ = ['aREA_classic']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def norm_ppf(x):
    return ndtri(x)

def sigT(x, slope = 20, inflection = 0.5):
    return (1 - 1/(1 + np.exp(slope * (x - inflection))))

def rankdata_torch_ordinal(x, device=None, verbose=False):
    if verbose:
        print(f"[rankdata_torch] using device={device}")

    with torch.inference_mode():
        # Move the input matrix to torch
        xt = torch.as_tensor(x, device=device)

        # Rank each row using the double-argsort trick
        order = torch.argsort(xt, dim=1)
        ranks = torch.argsort(order, dim=1) + 1

    dtype = np.float32 if device=="mps" else np.float64
    # Return ranks on CPU as float64 to match the rest of the code unless MPS
    return ranks.detach().cpu().numpy().astype(dtype)

def tail_prep_and_ndtri_torch(rankMat, gesInds, device=None, verbose=False):

    if verbose:
        print(f"[tail_prep_torch] using device={device}")

    with torch.inference_mode():
        dtype = torch.float32 if device == "mps" else torch.float64

        # Move the rank matrix and selected gene indices to torch
        rank_t = torch.as_tensor(rankMat, device=device, dtype=dtype)
        gesInds_t = torch.as_tensor(gesInds, device=device, dtype=torch.long)

        # Build the 2-tailed and 1-tailed matrices
        ges2T = rank_t / (rank_t.shape[1] + 1)
        ges1T = torch.abs(ges2T - 0.5) * 2
        ges1T = ges1T + (1 - torch.max(ges1T)) / 2

        # Keep only the intersecting genes
        ges2T = ges2T[:, gesInds_t]
        ges1T = ges1T[:, gesInds_t]

        # Clamp values away from 0 and 1 so icdf does not produce infinities
        eps = torch.finfo(dtype).eps
        ges2T = torch.clamp(ges2T, eps, 1 - eps)
        ges1T = torch.clamp(ges1T, eps, 1 - eps)


        # Handle the ICDF calculation (The problematic part for MPS)
        if device == "mps":
            # Move to CPU for icdf calculation
            ges2T_cpu = ges2T.cpu()
            ges1T_cpu = ges1T.cpu()

            normal = torch.distributions.Normal(
                torch.tensor(0.0, dtype=dtype), # Defaults to CPU
                torch.tensor(1.0, dtype=dtype)
            )

            ges2TQ = normal.icdf(ges2T_cpu).to(device)
            ges1TQ = normal.icdf(ges1T_cpu).to(device)
        else:
            # Convert probabilities to z-scores with the inverse normal CDF
            normal = torch.distributions.Normal(
                torch.tensor(0.0, device=device, dtype=dtype),
                torch.tensor(1.0, device=device, dtype=dtype),
            )
            ges2TQ = normal.icdf(ges2T)
            ges1TQ = normal.icdf(ges1T)

    return ges2TQ, ges1TQ

def enrichment_dots_torch_same_output(
    ic_mat,
    mor_mat,
    ges2TQ,
    ges1TQ,
    samplesIndex,
    device=None,
    verbose=False,
):

    if verbose:
        print(f"[enrichment_dots_torch] using device={device}")

    with torch.inference_mode():
        cpu_type = np.float32 if device == "mps" else np.float64
        torch_type = torch.float32 if device == "mps" else torch.float64

        # Move interaction confidence and mode of regulation to a tensor
        ic_t = torch.as_tensor(
            ic_mat.to_numpy(dtype=cpu_type, copy=False),
            device=device,
            dtype=torch_type,
        )
        mor_t = torch.as_tensor(
            mor_mat.to_numpy(dtype=cpu_type, copy=False),
            device=device,
            dtype=torch_type,
        )

        # Convert the directional and undirectional z scored to tensor
        ges2TQ = torch.as_tensor(ges2TQ, device=device, dtype=torch_type)
        ges1TQ = torch.as_tensor(ges1TQ, device=device, dtype=torch_type)

        # Directed enrichment: (IC * MoR)^T @ ges2TQ^T
        # Basically dES = sum ocer target genes if ic x mor x directional gene scores
        d_weights = (ic_t * mor_t).transpose(0, 1)
        d_scores = torch.matmul(d_weights, ges2TQ.transpose(0, 1))

        # Undirected enrichment: ((1 - abs(MoR)) * IC)^T @ ges1TQ^T
        # Basically uES = sum over target genes if ic x non-directional weight x extremeness
        u_weights = ((1.0 - torch.abs(mor_t)) * ic_t).transpose(0, 1)
        u_scores = torch.matmul(u_weights, ges1TQ.transpose(0, 1))

    # Convert back to pandas with the same labels as before
    dES = pd.DataFrame(
        d_scores.detach().cpu().numpy(),
        index=ic_mat.columns,
        columns=samplesIndex,
    )
    uES = pd.DataFrame(
        u_scores.detach().cpu().numpy(),
        index=ic_mat.columns,
        columns=samplesIndex,
    )

    return dES, uES

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def aREA_classic(gex_data, interactome, layer = None, eset_filter = False, min_targets=30, device="cpu", rank_ordinal=False, verbose=False, verbose_timing=False):
    """\
    Allows the individual to infer normalized enrichment scores from gene
    expression data using the analytical ranked enrichment analysis (aREA)
    function.

    It is the basis of the VIPER (Virtual Inference of Protein-activity
    by Enriched Regulon analysis) algorithm.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy) or in a
        pd.DataFrame.
    interactome
        The interactome object.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
    """

    if(verbose): print(f"Using {device} for aREA calculations")

    # If device is none, check if cuda is available. Otherwise, use cpu
    if device is None:
        if torch.cuda.is_available():
            device = "cuda"
        else:
            device = "cpu"

    # Set the timing variables/functions
    timings = {}
    t_total = time.perf_counter()

    # Helper function for marking the time
    def mark(name, t0):
        timings[name] = time.perf_counter() - t0

    # Filter out those with target less than min.targets
    interactome = interactome.copy()
    t0 = time.perf_counter()

    # Convert gex to dataframe
    if isinstance(gex_data, AnnData):
        gex_df = gex_data.to_df(layer)
        gex_data = None #deassign refrence
    elif isinstance(gex_data, pd.DataFrame):
        gex_df = gex_data
    else:
        raise ValueError("gex_data is type: " + str(type(gex_data)) + ". Must be anndata.AnnData or pd.DataFrame.")

    if (eset_filter):
        # This will affect the rankings of genes by eliminating those not present in the interactome
        tmp = np.unique(np.concatenate((
            interactome.get_target_names().astype(str),
            interactome.get_reg_names().astype(str)
        )))
        gex_df = gex_df.iloc[:,gex_df.columns.isin(pd.Series(tmp))]

    # Collect important data from gex_df before creating ges_arr
    varNames = gex_df.columns.to_list()
    samplesIndex = gex_df.index
    ges_arr = gex_df.values

    # Prune the interactome
    interactome.prune(min_targets = min_targets, max_targets = None, eliminate = False, verbose = False)

    # Mark the time for input prep and reset the counter
    mark("input_prep", t0)
    t0 = time.perf_counter()

    # ------------ find intersecting genes ------------
    # Use get_target_names as part of the interactome class to get a list of all targets in the interactome
    targetSet = interactome.get_target_names()

    # Convert to a set for fast lookup
    varNames_set = set(varNames)
    varname_to_idx = {g: i for i, g in enumerate(varNames)}

    # Get the intersction of gene names in the gExpr signature and those in the target set
    intersectGenes = np.array([x for x in targetSet if x in varNames_set])

    n_targets_not_in_exp_genes = np.count_nonzero(~np.isin(targetSet, varNames))
    if n_targets_not_in_exp_genes > 0:
        warnings.warn('interactome "' + str(interactome.name) + '" contains ' +
                        str(n_targets_not_in_exp_genes) + " targets missing from gex_data.var.\n\t" +
                        "Please run interactome.filter_targets(gex_data.var_names) on your network to\n\t" +
                        "resolve this. It is highly recommend to do this on the unPruned network and\n\t"+
                        "then prune. This way the Pruned network contains a consistent number of targets per\n\t"
                        "regulator, all of which exist within gex_data.")
        interactome.filter_targets(varNames, verbose = verbose)

    # Precompute gene indices here too
    gesInds = [varname_to_idx[i] for i in intersectGenes]

    # mark the time for filtering and reset the counter
    mark("intersect_and_filter_targets", t0)
    t0 = time.perf_counter()

    # rank transform the GES using either the rankdata function from scipy.stats for averaged ranking
    # or pytorch ranking for ordinal, which can be less stable/consistent but faster.
    # Using the ranks, it then assigns each gene a score based off of the inverse CDF for a
    # standard distribution (z-like score), so some genes can receive different value. The sign
    # of the NES is based soley off of the sign of the dES. Therefore, if the dES was already close
    # to 0, this small difference can have the effect of flipping the sign of some protein NES scores.
    if rank_ordinal:
        if(verbose): print("Rank transforming the data with ordinal ranking")
        rankMat = rankdata_torch_ordinal(ges_arr, device=device, verbose=verbose)
    else:
        if(verbose): print("Rank transforming the data with averaged ranking")
        rankMat = rankdata(ges_arr, axis = 1)

    # mark the time for ranking and reset the counter
    mark("rankdata", t0)
    t0 = time.perf_counter()

    # ------------ prepare the 1-tailed / 2-tailed matrices ------------
    # gesInds is a series of indices - the index of every target in the gExpr signature matrix
    # for each of the intersecting genes.
    # To get the one tailed and two tailed matrices, we normalize our rank values between 0 and 1
    # across each sample, thereby scaling our data across each sample. To do this, we divide each
    # row of the rankMat by the number of genes plus 1 to get the 2-tailed matrix.
    # For a one tailed test:
    #     (1) since each sample has a range of 0 to 1, we recenter our values at 0 for each
    #         sample by subtracting 0.5
    #     (2) take the absolute value so values are classified by their extremeness and not their sign
    #     (3) return the range of data from 0 to 0.5 to 0 to 1 by multiplying by 2
    # If the maximum value is less than 1, then we add half the difference between 1 and the max
    # value in the matrix to keep all values in the desired range of 0 to 1. This is an extra
    # normalization procedure. Finally, we apply the inverse CDF of the standard normal distribution
    # to obtain the z-like score matrices ges2TQ and ges1TQ. These steps arehandled by the
    # helper function tail_prep_and_ndtri_torch().

    if(verbose): print("Preparing the 1-tailed / 2-tailed matrices")

    ges2TQ, ges1TQ = tail_prep_and_ndtri_torch(
        rankMat,
        gesInds,
        device=device,
        verbose=verbose
    )
    del rankMat

    # mark the time for ndtri and reset the counter
    mark("tail_prep_and_ndtri", t0)
    t0 = time.perf_counter()

    # ------------ reduce regulon matrices ------------
    # The ic_mat is the matrix with regulators in the columns, targets in the rows and likelihood (weights) as values
    # (we filter to intersectGenes as targets by using .loc[intersectGenes])
    if(verbose): print("Computing the likelihood matrix")
    ic_mat = interactome.ic_mat()#.loc[intersectGenes]

    # mark the time for ic_mat and reset the counter
    mark("ic_mat", t0)
    t0 = time.perf_counter()

    # The morDict is the matrix with regulators in the columns, targets in the rows and tfmode (modes) as values
    if(verbose): print("Computing the modes matrix")
    mor_mat = interactome.mor_mat()#.loc[intersectGenes]

    # mark the time for mor_mat and reset the counter
    mark("mor_mat", t0)
    t0 = time.perf_counter()

    # ------------ directed / undirected enrichment ------------
    # For directed enrichment, we multiply the interaction confidence matrix (ic_mat)
    # and the tfmode matrix (mor_mat) to get directional weights of the targets in
    # each regulon. We then perform a dot product of these weights with the genes in
    # our 2-tailed z-score matrix ges2TQ to obtain the directed enrichment scores dES.

    # For undirected enrichment, we multiply the interaction confidence matrix by
    # (1 - abs(mor_mat)) to get non-directional weights of the targets in each regulon.
    # We then perform a dot product of these weights with the genes in our 1-tailed
    # z-score matrix ges1TQ to obtain the undirected enrichment scores uES.

    # These same calculations are now handled in PyTorch by the helper function
    # enrichment_dots_torch_same_output(), which returns the results in the same
    # pandas DataFrame format as before.

    if(verbose): print("Computing directed and undirected enrichment")

    dES, uES = enrichment_dots_torch_same_output(
        ic_mat=ic_mat,
        mor_mat=mor_mat,
        ges2TQ=ges2TQ,
        ges1TQ=ges1TQ,
        samplesIndex=samplesIndex,
        device=device,
        verbose=verbose,
    )
    mark("gpu_enrichment_dots_total", t0)

    del ges2TQ
    del ges1TQ
    del mor_mat
    del ic_mat

    t0 = time.perf_counter()

    # ------------ Integrate enrichment ------------
    if(verbose): print("Integrating enrichment")
    # We integrate our directed and undirected enrichment scores matrices to get our integrated enrichment scores
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)
    del dES
    del uES

    # ------------ make NES (Normalized Enrichment Scores) matrix ------------
    # interactome.icp_vec() returns a vector
    # This vector is generated by taking each individual regulon in the newtork and calculating
    # the likelihood index proportion to all interactions
    # The icProportion function that makes this calculation is in pyviper_classes.py.
    # This icProportion function calculates the proportion of the "Interaction Confidence" (IC) score
    # for each interaction in a network, relative to the maximum IC score in the network.
    # We multiply the Interaction Confidence of each regulator across the scores of each regulator
    # by making these multplications along the index (axis = 0)
    nES = iES.mul(interactome.icp_vec(), 0)
    del iES
    # We transpose our NES matrix to have samples in the rows and regulators in the columns
    nES = np.transpose(nES)
    # We make the rownames for the samples the same as that in the gExpr matrix
    nES.index = samplesIndex
    # Remove 'regulator' from the columns' name
    nES.columns.name = None

    # Mark the finalize nes time
    mark("finalize_nes", t0)

    # Print out the timings if the argument is true
    timings["total"] = time.perf_counter() - t_total

    # If verbose timing, print out all the time steps for each section
    if verbose_timing:
        print("\n[aREA timing]")
        total = timings["total"]
        for k, v in timings.items():
            pct = 100 * v / total if total > 0 else 0
            print(f"{k:>28s}: {v:8.4f} s   ({pct:5.1f}%)")

    # We return our result
    return(nES)
