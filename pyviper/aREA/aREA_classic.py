### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from scipy.stats import norm

### ---------- EXPORT LIST ----------
__all__ = ['aREA_classic']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def sigT(x, slope = 20, inflection = 0.5):
    return (1 - 1/(1 + np.exp(slope * (x - inflection))))

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def aREA_classic(gex_data, interactome, layer = None, eset_filter = False, min_targets=30, verbose = True):
    """\
    Allows the individual to infer normalized enrichment scores from gene
    expression data using the analytical ranked enrichment analysis (aREA)
    function.

    It is the basis of the VIPER (Virtual Inference of Protein-activity
    by Enriched Regulon analysis) algorithm.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy).
    interactome
        The interactome object.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
    """
    # Filter out those with target less than min.targets
    interactome = interactome.copy()
    interactome.prune(min_targets = min_targets, max_targets = None, eliminate = False, verbose = False)

    if (eset_filter):
        # This will affect the rankings of genes by eliminating those not present in the interactome
        tmp = np.unique(np.concatenate((interactome.get_targetSet(), interactome.get_regulonNames())))
        gex_data = gex_data[:,gex_data.var_names.isin(pd.Series(tmp))]

    if layer is None:
        gesMat = gex_data.X
    else:
        gesMat = gex_data.layers[layer]

    # ------------ find intersecting genes ------------
    # Use get_targetSet as part of the interactome class to get a list of all targets in the interactome
    targetSet = interactome.get_targetSet()
    # Get a list of the gene names of the gExpr signature matrix
    varNames = gex_data.var_names.to_list()
    # Get the intersction of gene names in the gExpr signature and those in the target set
    intersectGenes = np.intersect1d(targetSet, varNames)

    n_targets_not_in_exp_genes = np.count_nonzero(~np.isin(targetSet, varNames))
    if n_targets_not_in_exp_genes > 0:
        raise ValueError('interactome "' + str(interactome.name) + '" contains ' +
                         str(n_targets_not_in_exp_genes) + " targets missing from gex_data.var.\n\t" +
                        "Please run interactome.filter_targets(gex_data.var_names) on your network to\n\t" +
                         "resolve this. It is highly recommend to do this on the unPruned network and\n\t"+
                         "then prune to the pruned network contains a consistent number of targets per\n\t"
                         "regulator, allow of which exist within gex_data.")

    # rank transform the GES using the rankdata function from scipy.stats
    if(verbose): print("Rank transforming the data")
    rankMat = rankdata(gesMat, axis = 1)

    # ------------ reduce regulon matrices ------------
    # The icMat is the matrix with regulators in the columns, targets in the rows and likelihood (weights) as values
        # (we filter to intersectGenes as targets by using .loc[intersectGenes])
    if(verbose): print("Computing the likelihood matrix")
    icMat = interactome.icMat()#.loc[intersectGenes]
    # The morDict is the matrix with regulators in the columns, targets in the rows and tfmode (modes) as values
    if(verbose): print("Computing the modes matrix")
    morMat = interactome.morMat()#.loc[intersectGenes]

    # ------------ prepare the 1-tailed / 2-tailed matrices ------------
    if(verbose): print("Preparing the 1-tailed / 2-tailed matrices")
    # gesInds is a series of indices - the index of every target in the gExpr signature matrix
        # for each of the intersecting genes
    gesInds = [varNames.index(i) for i in intersectGenes]

    # To get the one tailed matrix, we normalize our rank values between 0 and 1 across each sample,
        # thereby scaling our data across each sample
    # To do this, we divide each row of the rankMat by the number of samples (rows) plus 1 to get the 2-tailed matrix
    ges2T = rankMat / (rankMat.shape[1] + 1)
    # For a one tailed test,
        # (1) since each sample has a range of 0 to 1, we recenter our values at 0 for each sample by subtracting 0.5
        # (2) take the absolute value so values are classified by their extremeness and not their sign
        # (3) return the range of data from 0 to 0.5 to 0 to 1 by multiplying by 2
    ges1T = abs(ges2T - 0.5) * 2
    # If the maximum value is less than 1, then we add half the difference between 1 and the max value in the matrix
        # to keep all values in the the desired range of 0 to 1. This is an extra normalization procedure
    ges1T = ges1T + (1 - np.max(ges1T))/2
    # norm.ppf is a function from the NumPy library that calculates the percent point function of the normal distribution,
    # which is the inverse of the cumulative distribution function. In this code, norm.ppf(ges2T[:, gesInds])
    # and norm.ppf(ges1T[:, gesInds]) calculate the z-scores for each column in gesInds, using the corresponding
    # values from ges2T and ges1T. The resulting matrices are called ges2TQ and ges1TQ, respectively.
    ges2TQ = norm.ppf(ges2T[:, gesInds])
    ges1TQ = norm.ppf(ges1T[:, gesInds])

    if(verbose): print("Computing enrichment")
    # ------------ 2-tail enrichment ------------
    # We multiply the likelihood matrix (icMat) and the tfmode matrix (morMat)
    # to get directional weights of the targets in each regulon
    dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
    # We then perform a dot product of these weights and the genes in our 2 tailed Z-scores matrix ges2TQ
    # to get our directed enrichment scores (samples in columns, regulators in the rows)
    dES = dES.dot(np.transpose(ges2TQ))

    # ------------ 1-tail enrichment ------------
    # We multiply the likelihood matrix (icMat) and the tfmode matrix (morMat)
    # to get undirected weights of the targets in each regulon
        # The farther the tfmode is from 0 and closer it is to 1, the smaller the weights
    uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
    # We then perform a dot product of these weights and the genes in our 1 tailed Z-scores matrix ges1TQ
    # to get our undirected enrichment scores (samples in columns, regulators in the rows)
    uES = uES.dot(np.transpose(ges1TQ))

    # ------------ Integrate enrichment ------------
    if(verbose): print("Integrating enrichment")
    # We integrate our directed and undirected enrichment scores matrices to get our integrated enrichment scores
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)

    # ------------ make NES (Normalized Enrichment Scores) matrix ------------
    # interactome.icpVec() returns a vector
        # This vector is generated by taking each individual regulon in the newtork and calculating
        # the likelihood index proportion to all interactions
        # The icProportion function that makes this calculation is in pyviper_classes.py.
            # This icProportion function calculates the proportion of the "Interaction Confidence" (IC) score
            # for each interaction in a network, relative to the maximum IC score in the network.
    # We multiply the Interaction Confidence of each regulator across the scores of each regulator
        # by making these multplications along the index (axis = 0)
    nES = iES.mul(interactome.icpVec(), 0)
    # We transpose our NES matrix to have samples in the rows and regulators in the columns
    nES = np.transpose(nES)
    # We make the rownames for the samples the same as that in the gExpr matrix
    nES.index = gex_data.obs.index
    # We return our result
    return(nES)

# def bootstrap_aREA(gesObj, intObj, bmean, bsd, eset_filter = False):
#     if (eset_filter):
#         tmp = list(set(list(intObj.get_targetSet())+list(intObj.get_regulonNames())))
#         gesObj = gesObj[:,gesObj.var_names.isin(pd.Series(tmp))]
#
#     samples = gesObj.shape[0]
#     regulons = len(intObj.get_regulonNames())
#
#     nes = np.zeros((samples,regulons))
#     nessd = np.zeros((samples,regulons))
#
#     bootstrap_obj = anndata.AnnData(X = np.zeros(bmean.shape),var=gesObj.var)
#
#     for i in range(samples): # for each row
#
#         bootstrap_obj.X = (-gesObj.X[i,:] + bmean)/bsd
#
#         result = aREA(bootstrap_obj, intObj, eset_filter = False)
#         nes[i] = np.mean(result, axis = 0)
#         nessd[i] = np.std(result, axis = 0)
#
#     result = pd.DataFrame(nes,index=gesObj.obs.index,columns=result.columns)
#
#     result = result.reset_index().melt(
#         id_vars = 'index',
#         var_name = 'gene')
#
#     return result
