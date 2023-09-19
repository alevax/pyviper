import numpy as np
from scipy.stats import norm
from scipy.stats import rankdata
from scipy.stats import ttest_1samp
import os
import shutil
import pandas as pd
import anndata
import scanpy as sc
from pyther_classes import *
from pyther_narnea import *
import pathlib
from tqdm import tqdm
from joblib import Parallel, delayed

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- CORE FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# calculator

def sigT(x, slope = 20, inflection = 0.5):
    return (1 - 1/(1 + np.exp(slope * (x - inflection))))

def sample_ttest(i,array):
    return ttest_1samp((array[i] - np.delete(array, i, 0)), 0).statistic



def load_interactome_from_tsv(filePath, intName):
    """\
    Allows the user to load an interactome object from a TSV file.

    The TSV file is created by the R function InteractomeToTable.

    Parameters
    ----------
    filePath
        The file to the regulon.tsv file.
    intName
        The name of the interactome.
    Returns
    -------
    A pd.DataFrame of :class:`~pyther_classes.Interactome`.
    """
    # read file
    netTable = pd.read_csv(filePath, sep = '\t')
    interactome = Interactome('intName', netTable) # what's the opint of this 'intname'
    # return
    return(interactome)

# def load_interactome_from_tsv(filePath, intName):
#     """\
#     Allows the user to load an interactome object from a TSV file.
#
#     The TSV file is created by the R function InteractomeToTable.
#
#     Parameters
#     ----------
#     filePath
#         The file to the regulon.tsv file.
#     intName
#         The name of the interactome.
#     Returns
#     -------
#     A dictionary of :class:`~pyther_classes.Interactome`.
#     """
#     # read file
#     netTable = pd.read_csv(filePath, sep = '\t')
#     interactome = Interactome('intName')
#     # loop through regulators
#     uniqueRegs = netTable.regulator.unique()
#     for u in tqdm(uniqueRegs, desc="Processing regulators", unit="regulator"):
#         # subset dataframe
#         uDF = netTable[netTable.regulator == u]
#         # make dictionaries
#         icDict = dict(zip(uDF.target, uDF.likelihood))
#         morDict = dict(zip(uDF.target, uDF.mor))
#         # make regulon object
#         regObj = Regulon(u, icDict, morDict)
#         interactome.addReg(u, regObj)
#     # return
#     return(interactome)

# def aREA(gex_data, interactome, layer = None):
#     """\
#     Allows the individual to infer normalized enrichment scores from gene
#     expression data using the analytical ranked enrichment analysis (aREA)
#     function.
#
#     It is the basis of the VIPER (Virtual Inference of Protein-activity
#     by Enriched Regulon analysis) algorithm.
#
#     Parameters
#     ----------
#     gex_data
#         Gene expression stored in an anndata object (e.g. from Scanpy).
#     interactome
#         The interactome object.
#     layer
#         The layer in the anndata object to use as the gene expression input
#         (default = None).
#     Returns
#     -------
#     A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
#     """
#     if layer is None:
#         gesMat = gex_data.X
#     else:
#         gesMat = gex_data.layers[layer]
#
#     # rank transform the GES
#     rankMat = rankdata(gesMat, axis = 1)
#
#     # find intersecting genes
#     targetSet = interactome.get_targetSet()
#     varNames = gex_data.var_names.to_list()
#     intersectGenes = [value for value in targetSet if value in varNames]
#
#     # reduce regulon matrices
#     icMat = interactome.icMat().loc[intersectGenes]
#     morMat = interactome.morMat().loc[intersectGenes]
#
#     # prepare the 1-tailed / 2-tailed matrices
#     gesInds = [varNames.index(i) for i in intersectGenes]
#     ges2T = rankMat / (rankMat.shape[1] + 1)
#     ges1T = abs(ges2T - 0.5) * 2
#     ges1T = ges1T + (1 - np.max(ges1T))/2
#     ges2TQ = norm.ppf(ges2T[:, gesInds])
#     ges1TQ = norm.ppf(ges1T[:, gesInds])
#
#     # 2-tail enrichment
#     dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
#     dES = dES.dot(np.transpose(ges2TQ))
#
#     # 1-tail enrichemnt
#     uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
#     uES = uES.dot(np.transpose(ges1TQ))
#
#     # integrate
#     iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)
#     # make NES
#     nES = iES.mul(interactome.icpVec(), 0)
#     nES = np.transpose(nES)
#     nES.index = gex_data.obs.index
#
#     return(nES)

def aREA(gex_data, interactome, eset_filter = False, layer = None):
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
    if (eset_filter):
        tmp = list(set(list(interactome.get_targetSet())+list(interactome.get_regulonNames())))
        gex_data = gex_data[:,gex_data.var_names.isin(pd.Series(tmp))]
        #would this line have influnence outside this function?

    if layer is None:
        gesMat = gex_data.X
    else:
        gesMat = gex_data.layers[layer]



    # rank transform the GES using the rankdata function from scipy.stats
    rankMat = rankdata(gesMat, axis = 1)

    # ------------ find intersecting genes ------------
    # Use get_targetSet as part of the interactome class to get a list of all targets in the interactome
    targetSet = interactome.get_targetSet()
    # Get a list of the gene names of the gExpr signature matrix
    varNames = gex_data.var_names.to_list()
    # Get the intersction of gene names in the gExpr signature and those in the target set
    intersectGenes = [value for value in targetSet if value in varNames]

    # ------------ reduce regulon matrices ------------
    # The icMat is the matrix with regulators in the columns, targets in the rows and likelihood (weights) as values
        # (we filter to intersectGenes as targets by using .loc[intersectGenes])
    icMat = interactome.icMat().loc[intersectGenes]
    # The morDict is the matrix with regulators in the columns, targets in the rows and tfmode (modes) as values
    morMat = interactome.morMat().loc[intersectGenes]

    # ------------ prepare the 1-tailed / 2-tailed matrices ------------

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

    # ------------ 2-tail enrichment ------------
    # We multiply the likelihood matrix (icMat) and the tfmode matrix (morMat)
    # to get directional weights of the targets in each regulon
    dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
    # We then perform a dot product of these weights and the genes in our 2 tailed Z-scores matrix ges2TQ
    # to get our directed enrichment scores (samples in columns, regulators in the rows)
    dES = dES.dot(np.transpose(ges2TQ))

    # ------------ 1-tail enrichemnt ------------
    # We multiply the likelihood matrix (icMat) and the tfmode matrix (morMat)
    # to get undirected weights of the targets in each regulon
        # The farther the tfmode is from 0 and closer it is to 1, the smaller the weights
    uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
    # We then perform a dot product of these weights and the genes in our 1 tailed Z-scores matrix ges1TQ
    # to get our undirected enrichment scores (samples in columns, regulators in the rows)
    uES = uES.dot(np.transpose(ges1TQ))

    # ------------ Integrate enrichment ------------
    # We integrate our directed and undirected enrichment scores matrices to get our integrated enrichment scores
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)

    # ------------ make NES (Normalized Enrichment Scores) matrix ------------
    # interactome.icpVec() returns a vector
        # This vector is generated by taking each individual regulon in the newtork and calculating
        # the likelihood index proportion to all interactions
        # The icProportion function that makes this calculation is in pyther_classes.py.
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

def bootstrap_aREA(gesObj, intObj, bmean, bsd, eset_filter = False):
    if (eset_filter):
        tmp = list(set(list(intObj.get_targetSet())+list(intObj.get_regulonNames())))
        gesObj = gesObj[:,gesObj.var_names.isin(pd.Series(tmp))]

    samples = gesObj.shape[0]
    regulons = len(intObj.get_regulonNames())

    nes = np.zeros((samples,regulons))
    nessd = np.zeros((samples,regulons))

    bootstrap_obj = anndata.AnnData(X = np.zeros(bmean.shape),var=gesObj.var)

    for i in range(samples): # for each row

        bootstrap_obj.X = (-gesObj.X[i,:] + bmean)/bsd

        result = aREA(bootstrap_obj, intObj, eset_filter = False)
        nes[i] = np.mean(result, axis = 0)
        nessd[i] = np.std(result, axis = 0)

    result = pd.DataFrame(nes,index=gesObj.obs.index,columns=result.columns)

    result = result.reset_index().melt(
        id_vars = 'index',
        var_name = 'gene')

    return result

def meta_aREA(gesObj, intObj, eset_filter = False, pleiotropy = False, pleiotropyArgs = {}, layer = None, mvws = 1, njobs = 1, verbose = False):
    if type(intObj) == Interactome:
        preOp = aREA(gesObj, intObj, eset_filter, layer)
    
    #for I've added 'preparing networks' after input into main function.
    #if len(intObj) == 1:
    #    preOp = aREA(gesObj, intObj[0], eset_filter, layer)

    elif njobs == 1:
        netMets = [aREA_melt(gesObj, iObj, eset_filter, pleiotropy, pleiotropyArgs, layer) for iObj in intObj]
        preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)
    else:
        joblib_verbose = 0
        if verbose:
            print("Computing regulons enrichment with aREA")
            joblib_verbose = 11
        # n_jobs need to be decided.
        netMets = Parallel(n_jobs = njobs, verbose = joblib_verbose)(
            (delayed)(aREA_melt)(gesObj, iObj, eset_filter, pleiotropy, pleiotropyArgs, layer)
            for iObj in intObj
            )
        preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)
    return preOp

def aREA_melt(gesObj, intObj, eset_filter = False, pleiotropy = False, pleiotropyArgs = {}, layer = None,):
    pb = None
    result = aREA(gesObj, intObj, eset_filter, layer)
    result = result.reset_index().melt(
        id_vars = 'index',
        var_name = 'gene'
    )

    if (pleiotropy):
        print('pleiotropy is currently unavailable')

    return result

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ PYTHER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
def mat_to_anndata(mat):
    # Helper function for *pyther* and *path_enr*
    # Create obs dataframe
    mat_sampleNames = pd.DataFrame(index=range(len(mat.index.values)),columns=range(0))
    mat_sampleNames.index = mat.index.values
    mat_sampleNames

    # Create var dataframe
    mat_features = pd.DataFrame(index=range(len(mat.columns.values)),columns=range(0))
    mat_features.index = mat.columns.values
    mat_features

    # Convert the pandas dataframe from Pyther into a new Anndata object
    pax_data = anndata.AnnData(X=mat,
                               obs=mat_sampleNames,
                               var=mat_features)
    return(pax_data)

def consolidate_meta_aREA_results(netMets, mvws = 1, verbose = False):
    firstMat = netMets.pop(0)

    for thisMat in netMets:
        firstMat = firstMat.merge(thisMat,how = 'outer',on = ['index','gene'])

    firstMat.fillna(0,inplace = True)

    result = firstMat[['index','gene']]
    nes = firstMat[list(firstMat.columns)[2:]].values

    #mvws = 1
    if type(mvws) == int :
        ws = np.abs(nes)**mvws
        if verbose:
            print('mvws =' , mvws)
    else:
        ws = sigT(np.abs(nes),mvws[1],mvws[0])

    result['value'] = np.sum(nes*ws,axis =1)/np.sum(ws,axis =1)

    preOp = result.pivot(index='index',columns="gene", values="value")
    return preOp

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------
def pyther(gex_data,
           interactome,
           layer = None,
           njobs = 1,
           eset_filter = True,
           # bootstrap = 0,
           # dnull = None,
           # pleiotropy = False,
           # minsize=25,
           # adaptive_size=False,
           method=[None, "scale", "rank", "mad", "ttest"],
           enrichment = [None, 'area','narnea'],
           mvws=1,
           # pleiotropyArgs={'regulator':0.05, 'shadow':0.05, 'targets':10, "penalty":20, "method":"adaptive"},
           verbose = True,
           output_type  = ['anndata', 'ndarray']
           ):

#
    gesObj = gex_data
    intList = interactome
    pd.options.mode.chained_assignment = None

    # if verbose:
    #     print('preparing regulon networks')

    # if type(interactome) == pd.core.frame.DataFrame:
    #     intList = Interactome('intName', interactome)

    # else:

    #     for i in interactome:
    #         intList.append(Interactome('intName', i))


    #adapt for single regulon input is nolonger needed.

    # if type(intList) == Interactome:
    #     intList = [intList]
    #     print('Single regulon inputed')

    # if pleiotropy:
    #     bootstrap = 0
    #     if verbose:
    #         print("Using pleiotropic correction, bootstraps iterations are ignored.")

    if verbose:
        print("Preparing the association scores")

    ''' not in current version 3.8 of python , only in python 3.10
    match method:

        case 'scale': #scale each gene in all sample
            gesObj.X = (gesObj.X - np.mean(gesObj.X,axis=0))/np.std(gesObj.X,axis=0)
        case 'rank':
            gesObj.X = rankdata(gesObj.X,axis=0)*(np.random.random(gesObj.X.shape)*2/10-0.1)
        case 'mad':
            median = np.median(gesObj.X,axis=0)
            gesObj.X = (gesObj.X-median)/np.median(np.abs(gesObj.X-median))
        case 'ttest':
            # I dont quite understand this part maybe ask later:
            gesObj.X = gesObj.X
    '''
    if method == 'scale':
        gesObj.X = (gesObj.X - np.mean(gesObj.X,axis=0))/np.std(gesObj.X,axis=0)
    elif method == 'rank':
        gesObj.X = rankdata(gesObj.X,axis=0)*(np.random.random(gesObj.X.shape)*2/10-0.1)
    elif method == 'mad':
        median = np.median(gesObj.X,axis=0)
        gesObj.X = (gesObj.X-median)/np.median(np.abs(gesObj.X-median))
    elif method == 'ttest':
        gesObj.X = np.array([sample_ttest(i, gesObj.X.copy()) for i in range(gesObj.shape[0])])

#     if bootstrap > 0:
#
# #        num_sample = gesObj.X.shape[0]
#         bmean = np.zeros((bootstrap,gesObj.X.shape[1]))
#         bsd = np.zeros((bootstrap,gesObj.X.shape[1]))
#
#
#         sampled_indices = np.random.choice(gesObj.obs.index,bootstrap*len(gesObj.obs))
#
#         for i in range(bootstrap):
#             sample_name = sampled_indices[i*len(gesObj.obs):(i+1)*len(gesObj.obs)]
#             #sample mean
#             bmean[i] = np.mean(gesObj[sample_name,:].X, axis=0)
#             bsd[i] = np.std(gesObj[sample_name,:].X, axis=0)
#
#         #targets = intObj.get_targetSet()
#         # may need a bootstrap area
#         netMets = Parallel(n_jobs = njobs)(
#         (delayed)(bootstrap_aREA)(gesObj,iObj,eset_filter = True,layer = layer)
#         for iObj in intList
#         )
#
#         preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)
#
#     else:



    if enrichment is None: enrichment = 'narnea'

    if enrichment == 'area':
        preOp = meta_aREA(gesObj, intList, eset_filter, layer = layer, njobs = njobs, verbose = verbose, mvws = mvws )
        if output_type == 'ndarray':
            op = preOp
        else:
            op = mat_to_anndata(preOp)
    else:
        # time record disabled
        #joblib_verbose = 0
        if verbose:
            print("Computing regulons enrichment with NaRnEa")
            #joblib_verbose = 11
        preOp = meta_narnea(gesObj, intList, sample_weight = True, njobs = njobs, verbose = verbose)
        if output_type == 'ndarray':
            op = preOp
        else:
            op = mat_to_anndata(preOp["nes"]) #anndata.AnnData(preOp["nes"])
            op.layers['pes'] = preOp["pes"]

    return op #final result

def path_enr(adata,
             interactome,
             layer = None,
             verbose = False):
    """\
    Allows the individual to infer normalized enrichment scores of pathways
    using the analytical ranked enrichment analysis (aREA) function.

    This is a wrapper for the aREA function in Pyther.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    interactome
        The interactome object containing pathways.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
    """
    if(verbose): print("Checking interactome names format...")
    interactome_gene_name_format = detect_interactome_name_type(interactome)
    if(verbose): print("Interactome names are formatted as " + interactome_gene_name_format + ".")
    if(verbose): print("Checking adata names format...")
    adata_gene_name_format_original = detect_index_name_type(adata)
    if(verbose): print("adata names are formatted as " + adata_gene_name_format_original + ".")
    make_adata_names_format_match_interactome = False

    if(interactome_gene_name_format != adata_gene_name_format_original):
        make_adata_names_format_match_interactome = True
        if(verbose): print("Translating adata names to match interactome...")
        adata = translate_adata_index(adata,
                                      current_format = adata_gene_name_format_original,
                                      desired_format = interactome_gene_name_format)

    # aREA takes the pathways interactome and the adata
    if(verbose): print("Running aREA using to calculate pathway enrichment...")
    path_enr_mat = aREA(adata, interactome, layer)

    if(make_adata_names_format_match_interactome is True):
        if(verbose): print("Returning adata names to original state...")
        adata = translate_adata_index(adata,
                                      current_format = interactome_gene_name_format,
                                      desired_format = adata_gene_name_format_original)
    # Create a new Anndata object
    pwe_data = mat_to_anndata(path_enr_mat)
    # This means we did pathway enrichment on VIPER: adata is pax_data
    if hasattr(adata, "gex_data"):
        pwe_data.gex_data = adata.gex_data
        adata.gex_data = None
        pwe_data.pax_data = adata
    # This means we did pathway enrichment on gex: adata is gex_data
    else:
        pwe_data.gex_data = adata
    return(pwe_data)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ----------------------------- TRANSLATE FUNCTIONS ---------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# ---------------------------- TRANSLATE LOAD DATA FUNCS ---------------------------
def load_transate_csv(path_to_csv):
    with open(path_to_csv) as temp_file:
        regulator_set = [line.rstrip('\n') for line in temp_file]
    return(regulator_set)
def load_mouse2human():
    mouse2human = pd.read_csv(get_pyther_dir() + "/data/translate/human2mouse.csv")
    del mouse2human[mouse2human.columns[0]]
    return(mouse2human)
def load_human2mouse():
    human2mouse = pd.read_csv(get_pyther_dir() + "/data/translate/mouse2human.csv")
    del human2mouse[human2mouse.columns[0]]
    return(human2mouse)
# ---------------------------- DETECT FORMAT TYPE FUNCS ---------------------------
def detect_interactome_name_type(interactome):
    return(detect_name_type(np.array(list(interactome.get_targetSet()))))
def detect_index_name_type(adata):
    return(detect_name_type(adata.var.index))
def detect_name_type(input_array):
    nrow = len(input_array)
    if nrow == 0: return(None) #redundant: case handled by 0-length for loop
    found_match = False
    human2mouse = load_human2mouse()
    gene_name_format = None
    i = 0
    for i in range(nrow):
        gene = input_array[i]
        if(gene in human2mouse["mouse_symbol"].values):
            gene_name_format = "mouse_symbol"
            break
        elif(gene in human2mouse["human_symbol"].values):
            gene_name_format = "human_symbol"
            break
        elif(gene in human2mouse["mouse_ensembl"].values):
            gene_name_format = "mouse_ensembl"
            break
        elif(gene in human2mouse["human_ensembl"].values):
            gene_name_format = "human_ensembl"
            break
    return(gene_name_format)
# ---------------------------- TRANSLATE SMALL PICTURE FUNCS ---------------------------
def translate_gene_name(input_gene_name,
              translate_df,
              current_format = "mouse_symbol",
              desired_format = "human_symbol"):
    translate_df_input_gene_row = translate_df.loc[translate_df[current_format] == input_gene_name]
    if(translate_df_input_gene_row.empty):
        output_gene_name = None
    else:
        output_gene_name = translate_df_input_gene_row[desired_format].values[0]
    return(output_gene_name)
def translate_adata_index_into_new_var_column(
         adata,
         translate_df,
         current_format = "mouse_symbol",
         desired_format = "human_symbol",
         show_progress_bar = True):
    current_gene_names = adata.var.index.values
    n_genes = len(current_gene_names)
    desired_gene_names = [None]*n_genes

    iterator = (i for i in range(n_genes)) if not show_progress_bar else tqdm(range(n_genes))
    for i in iterator:
        desired_gene_names[i] = translate_gene_name(
            current_gene_names[i],
            translate_df,
            current_format,
            desired_format
        )
    adata.var[desired_format] = desired_gene_names
    return(adata)
# ---------------------------- TRANSLATE BIG PICTURE FUNCS ---------------------------
def translate_adata_index_from_to_with_translate_df(adata,
         translate_df,
         current_format = "mouse_symbol",
         desired_format = "human_symbol",
         show_progress_bar = True):
    adata = translate_adata_index_into_new_var_column(
         adata,
         translate_df,
         current_format,
         desired_format,
         show_progress_bar)
    adata.var[current_format] = adata.var.index.values
    adata.var.set_index(desired_format, inplace=True)
    return(adata)
def translate_adata_index_from_to(adata,
                                  current_format = "mouse_symbol",
                                  desired_format = "human_symbol",
                                  show_progress_bar = True):
    acceptable_formats = ["mouse_symbol", "mouse_ensembl", "human_symbol", "human_ensembl"]
    if(current_format in ["mouse_symbol", "mouse_ensembl"]):
        translate_df = load_mouse2human()
    elif(current_format in ["human_symbol", "human_ensembl"]):
        translate_df = load_human2mouse()
    if current_format not in acceptable_formats:
        raise ValueError("Error: index of adata.var is not one the following:"
                         + "\n\t\t mouse_symbol, mouse_ensembl, human_symbol, human_ensembl")
    if(desired_format in adata.var.columns):
        adata.var[current_format] = adata.var.index
        adata.var.set_index(desired_format, inplace=True)
    else:
        adata = translate_adata_index_from_to_with_translate_df(adata,
             translate_df,
             current_format,
             desired_format = desired_format,
             show_progress_bar = show_progress_bar)
    return(adata)
def translate_adata_index(adata,
                          desired_format = "human_symbol",
                          current_format = None,
                          show_progress_bar = True):
    acceptable_formats = ["mouse_symbol", "mouse_ensembl", "human_symbol", "human_ensembl"]
    if desired_format not in acceptable_formats:
        raise ValueError("Error: desired_format is not one the following:"
                         + "\n\t\t mouse_symbol, mouse_ensembl, human_symbol, human_ensembl")
    if current_format is None:
        current_format = detect_index_name_type(adata)
    adata = translate_adata_index_from_to(adata,
                                          current_format,
                                          desired_format,
                                          show_progress_bar)
    return(adata)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# --------------------------- LOADING DATA FUNCTIONS --------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
def slice_concat(inner_function, gex_data ,bins = 10, write_local = True, **kwargs): 
    #kwargs are the parameters for the inner function.
    #slice the data cells * genes

    result_list = []
    size = int(gex_data.shape[0]/bins) 
    residue = gex_data.shape[0] % bins

    if write_local:
        os.mkdir('temp')

        for i in range(bins-1):
            segment = gex_data[i*size: i*size + size,]
            temp_result = inner_function(segment, **kwargs)

            if type(temp_result) == anndata._core.anndata.AnnData:
                temp_result = temp_result.to_df()


            temp_result.to_csv('temp/'+ str(i) + '.csv')
        
        # the last one
        segment = gex_data[(bins-1)*size: bins*size + residue,]
        temp_result = inner_function(segment, **kwargs)

        if type(temp_result) == anndata._core.anndata.AnnData:
            temp_result = temp_result.to_df()
        
        temp_result.to_csv('temp/'+ str(bins-1) + '.csv')


        all_file_list=os.listdir('temp')
        for single_file in all_file_list:
            result_list.append(pd.read_csv(os.path.join('temp',single_file)))

        shutil.rmtree('temp')

    else:        
        for i in range(bins):
            segment = gex_data[i*size: i*size + size,]
            result_list.append(inner_function(segment, **kwargs))

    
    # concat result

    result = pd.concat(result_list,axis=0).reset_index(drop = True)
    result.set_index(keys='index',inplace=True)
    return result

def get_pyther_dir():
    pyther_dir = str(pathlib.Path(__file__).parent.parent.resolve())
    return(pyther_dir)

def load_regulators(path_to_txt):
    with open(path_to_txt) as temp_file:
        regulator_set = [line.rstrip('\n') for line in temp_file]
    return(regulator_set)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

# ------------------------------- LOAD REGULATORS -----------------------------
def load_TFs(path_to_tfs = None, species = "human"):
    if path_to_tfs is None:
        if species == "human":
            path_to_tfs = get_pyther_dir() + "/data/regulatorIDs/tfs-hugo.txt"
        elif species == "mouse":
            path_to_tfs = get_pyther_dir() + "/data/regulatorIDs/tfs-hugo-mouse.txt"
    tfs_list = load_regulators(path_to_tfs)
    return(tfs_list)

def load_coTFs(path_to_cotfs = None, species = "human"):
    if path_to_cotfs is None:
        if species == "human":
            path_to_cotfs = get_pyther_dir() + "/data/regulatorIDs/cotfs-hugo.txt"
        elif species == "mouse":
            path_to_cotfs = get_pyther_dir() + "/data/regulatorIDs/cotfs-hugo-mouse.txt"
    cotfs_list = load_regulators(path_to_cotfs)
    return(cotfs_list)

def load_sig(path_to_sig = None, species = "human"):
    if path_to_sig is None:
        if species == "human":
            path_to_sig = get_pyther_dir() + "/data/regulatorIDs/sig-hugo.txt"
        elif species == "mouse":
            path_to_sig = get_pyther_dir() + "/data/regulatorIDs/sig-hugo-mouse.txt"
    sig_list = load_regulators(path_to_sig)
    return(sig_list)

def load_surf(path_to_surf = None, species = "human"):
    if path_to_surf is None:
        if species == "human":
            path_to_surf = get_pyther_dir() + "/data/regulatorIDs/surface-hugo.txt"
        elif species == "mouse":
            path_to_surf = get_pyther_dir() + "/data/regulatorIDs/surface-hugo-mouse.txt"
    surf_list = load_regulators(path_to_surf)
    return(surf_list)

# ---------------------------- LOAD MSIGDB REGULONS ---------------------------
def merge_dicts(dict1, dict2):
    res = {**dict1, **dict2}
    return res
def get_msigdb_reg_path(collection):
    reg_path = get_pyther_dir() + "/data/regulons/msigdb-" + collection + "-as-regulon.tsv"
    return(reg_path)
def load_msigdb_from_tsv(collection):
    reg_path = get_msigdb_reg_path(collection.lower())
    reg = load_interactome_from_tsv(reg_path, "MSigDB_" + collection.upper())
    return(reg)
def load_msigdb_regulon_single(collection = "c2"):
    reg = None
    if(collection.lower() in ["c2", "c5", "c6", "c7", "h"]):
        reg = load_msigdb_from_tsv(collection)
    return(reg)
def load_msigdb_regulon_multiple(collection = ["h", "c2"]):
    combined_dict = {}
    for i in range(len(collection)):
        new_dict = load_msigdb_regulon_single(collection[i]).regDict
        combined_dict = merge_dicts(combined_dict, new_dict)
    combined_interactome = Interactome(name = '+'.join(collection),
                                       regDict = combined_dict)
    return(combined_interactome)
def load_msigdb_regulon(collection = "h"):
    reg = None
    if(type(collection) is str):
        reg = load_msigdb_regulon_single(collection)
    elif(type(collection) is list):
        reg = load_msigdb_regulon_multiple(collection)
    return(reg)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ----------------------------- FILTERING FUNCTIONS ---------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------

# Python program to illustrate the intersection
# of two lists in most simple way
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def get_mat_from_anndata(adata, features_indices):
    mat = pd.DataFrame(adata.X[:,features_indices],
                       index = list(adata.obs_names),
                       columns = adata.var_names[features_indices,])
    return(mat)

def get_features_indices(adata, features_list):
    true_false_list = pd.Series(adata.var_names).isin(features_list).tolist()
    features_indices = [i for i, x in enumerate(true_false_list) if x]
    return(features_indices)

def get_features_list(adata,
                    feature_groups="all",
                    path_to_tfs = None,
                    path_to_cotfs = None,
                    path_to_sig = None,
                    path_to_surf = None):
    if(type(feature_groups) is str):
        feature_groups = [feature_groups]
    if "all" not in feature_groups:
        feature_groups = [x.lower() for x in feature_groups]
        features_list = list()
        if "tfs" in feature_groups or "tf" in feature_groups:
            tfs = load_TFs(path_to_tfs)
            features_list.extend(tfs)
        if "cotfs" in feature_groups or "cotf" in feature_groups:
            cotfs = load_coTFs(path_to_cotfs)
            features_list.extend(cotfs)
        if "sigs" in feature_groups or "sig" in feature_groups:
            sig = load_sig(path_to_sig)
            features_list.extend(sig)
        if "surfs" in feature_groups or "surf" in feature_groups:
            surf = load_surf(path_to_surf)
            features_list.extend(surf)
        features_list = intersection(features_list, list(adata.var_names))
    else:
        features_list = adata.var_names
    return(features_list)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

# ------------------------- GET ANNDATA OBJECT FILTERED -----------------------
def get_anndata_filtered_by_feature_group(adata,
                               feature_groups="all", #["TFs", "CoTFs", "sig", "surf"],
                               path_to_tfs = None,
                               path_to_cotfs = None,
                               path_to_sig = None,
                               path_to_surf = None):
    features_list = get_features_list(adata,
            feature_groups,
            path_to_tfs,
            path_to_cotfs,
            path_to_sig,
            path_to_surf)
    adata_filt = get_anndata_filtered_by_feature_list(adata, features_list)
    return(adata_filt)

def get_anndata_filtered_by_feature_list(adata, features_list):
    features_indices = get_features_indices(adata, features_list)
    mat = get_mat_from_anndata(adata, features_indices)
    adata_with_features_only = mat_to_anndata(mat)
    return(adata_with_features_only)

# ----------------------------- FILTER INTERACTOME ----------------------------
def prune_interactome(interactome, features_list):
    # First prune the regulators
    # Using a list comprehension to make a list of the keys to be deleted
    regulators_to_delete = [key for key in interactome.regDict if key not in features_list]
    # delete the key/s
    for regulator in regulators_to_delete:
        del interactome.regDict[regulator]

    return(interactome)

# ------------------------ SCANPY TOOLS PYTHER WRAPPERS -----------------------
def tl_pca(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]
    return(adata)

def tl_tsne(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.tsne(adata_filt, **kwargs)
    adata.obsm["X_tsne"] = adata_filt.obsm["X_tsne"]
    return(adata)

def tl_umap(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            svd_solver='arpack',
            n_neighbors=10,
            n_pcs=40,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    # Create a second anndata object
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # Move sc.pp.neighbors results to adata_filt
    # adata_filt.obsp["distances"] = adata.obsp["distances"]
    # adata_filt.obsp["connectivities"] = adata.obsp["connectivities"]
    # adata_filt.uns["neighbors"] = adata.uns["neighbors"]
    # Compute the UMAP
    sc.tl.umap(adata_filt, **kwargs)
    # Give the UMAP from the second anndata object to the original
    adata.obsm["X_umap"] = adata_filt.obsm["X_umap"]
    return(adata)

def tl_draw_graph(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.draw_graph(adata_filt, **kwargs)
    if "X_draw_graph_fa" in list(adata.obsm.keys()):
        adata.obsm["X_draw_graph_fa"] = adata_filt.obsm["X_draw_graph_fa"]
    if "X_draw_graph_fr" in list(adata.obsm.keys()):
        adata.obsm["X_draw_graph_fr"] = adata_filt.obsm["X_draw_graph_fr"]
    return(adata)

def tl_diffmap(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.diffmap(adata_filt, **kwargs)
    adata.obsm["X_diffmap"] = adata_filt.obsm["X_diffmap"]
    adata.uns['diffmap_evals'] = adata_filt.uns['diffmap_evals']
    return(adata)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ PLOTTING FUNCTIONS ---------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
def get_gex_anndata_with_nes_umap(adata):
    adata_gex = adata.gex_data
    adata_gex.obsm["X_umap"] = adata.obsm["X_umap"]
    return(adata_gex)

def get_pax_anndata_with_path_enr_umap(adata):
    adata_pax = adata.pax_data
    adata_pax.obsm["X_umap"] = adata.obsm["X_umap"]
    return(adata_pax)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------
def pl_umap(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.umap.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.umap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.umap(adata, **kwargs)

def pl_scatter(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.scatter.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.scatter function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.scatter(adata,**kwargs)

def pl_heatmap(adata,
     *,
     plot_stored_gex_data=False,
     plot_stored_pax_data=False,
     **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.heatmap.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.heatmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.heatmap(adata,**kwargs)

def pl_dotplot(adata,
                      *,
                      plot_stored_gex_data = False,
                      plot_stored_pax_data = False,
                      **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.dotplot(adata,**kwargs)

def pl_tracksplot(adata,
                      *,
                      plot_stored_gex_data = False,
                      plot_stored_pax_data = False,
                      **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.tracksplot(adata,**kwargs)

def pl_violin(adata,
                     *,
                     plot_stored_gex_data = False,
                     plot_stored_pax_data = False,
                     **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.violin(adata,**kwargs)

def pl_stacked_violin(adata,
                             *,
                             plot_stored_gex_data = False,
                             plot_stored_pax_data = False,
                             **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.stacked_violin(adata,**kwargs)

def pl_matrixplot(adata,
                         *,
                         plot_stored_gex_data = False,
                         plot_stored_pax_data = False,
                         **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.matrixplot(adata,**kwargs)

def pl_clustermap(adata,
                         *,
                         plot_stored_gex_data = False,
                         plot_stored_pax_data = False,
                         **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.clustermap(adata,**kwargs)

def pl_ranking(adata,
                      *,
                      plot_stored_gex_data = False,
                      plot_stored_pax_data = False,
                      **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.ranking(adata,**kwargs)

def pl_dendrogram(adata,
                         *,
                         plot_stored_gex_data = False,
                         plot_stored_pax_data = False,
                         **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gex_anndata_with_nes_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pax_anndata_with_path_enr_umap(adata)
    sc.pl.dendrogram(adata,**kwargs)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# --------------------------- DATA TRANSFORM FUNCTIONS ------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------

def mad_from_R(x, center=None, constant=1.4826, low=False, high=False):
    if center is None:
        center=np.median(x)
    x = x[~np.isnan(x)] if np.isnan(x).any() else x
    n = len(x)
    if (low or high) and n % 2 == 0:
        if low and high:
            raise ValueError("'low' and 'high' cannot be both True")
        n2 = n // 2 + int(high)
        return constant * np.sort(np.abs(x - center))[n2]
    return constant * np.median(np.abs(x - center))

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

# Function assumes features as rows and observations as columns
# Numerator Functions:
    # Median - numpy.median
    # Mean - numpy.mean
# Denominator Functions:
    # Median absolute deviation - mad_from_R
    # Standard deviation - statistics.stdev
def rank_norm(x, NUM_FUN=np.median, DEM_FUN = mad_from_R, trim=0):
    rank = rankdata(x, axis=0)
    median = NUM_FUN(rank, axis=1, keepdims=True)#np.median(rank, axis=1, keepdims=True)
    mad = np.apply_along_axis(DEM_FUN, 1, rank)

    x = ((rank - median)/mad[:, np.newaxis])

    print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
    print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
    print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))
    print("- Features to Remove:")

    # Take care of infinite values
    max_finite = np.nanmax(x[np.isfinite(x)])
    min_finite = np.nanmin(x[np.isfinite(x)])
    x[np.isposinf(x)] = max_finite
    x[np.isneginf(x)] = min_finite

    x = np.where(np.isnan(x), np.nanmin(x), x)
    x = np.clip(x, a_min=np.nanmin(x), a_max=np.nanmax(x))
    print("- Removing NULL/NA features ...")
    x = x[~np.isnan(x).any(axis=1)]

    print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
    print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
    print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))

    return x
