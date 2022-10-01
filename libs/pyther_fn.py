from pyther_classes import *
import numpy as np
from scipy.stats import norm
from scipy.stats import rankdata
import pandas as pd
import anndata

def InteractomefromTSV(filePath, intName):
    # read file
    netTable = pd.read_csv(filePath, sep = '\t')
    intObj = Interactome('intName')
    # loop through regulators
    uniqueRegs = netTable.regulator.unique()
    for u in uniqueRegs:
        # subset dataframe
        uDF = netTable[netTable.regulator == u]
        # make dictionaries
        icDict = dict(zip(uDF.target, uDF.likelihood))
        morDict = dict(zip(uDF.target, uDF.mor))
        # make regulon object
        regObj = Regulon(u, icDict, morDict)
        intObj.addReg(u, regObj)
    # return
    return(intObj)

def aREA(gesObj, intObj):
    gesMat = gesObj.X
        
    # rank transform the GES
    rankMat = rankdata(gesMat, axis = 1)
        
    # find intersecting genes
    targetSet = intObj.get_targetSet()
    varNames = gesObj.var_names.to_list()
    intersectGenes = [value for value in targetSet if value in varNames]
        
    # reduce regulon matrices
    icMat = intObj.icMat().loc[intersectGenes]
    morMat = intObj.morMat().loc[intersectGenes]
        
    # prepare the 1-tailed / 2-tailed matrices
    gesInds = [varNames.index(i) for i in intersectGenes]
    ges2T = rankMat / (rankMat.shape[1] + 1)
    ges1T = abs(ges2T - 0.5) * 2
    ges1T = ges1T + (1 - np.max(ges1T))/2
    ges2TQ = norm.ppf(ges2T[:, gesInds])
    ges1TQ = norm.ppf(ges1T[:, gesInds])
        
    # 2-tail enrichment
    dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
    dES = dES.dot(np.transpose(ges2TQ))
        
    # 1-tail enrichemnt
    uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
    uES = uES.dot(np.transpose(ges1TQ))
        
    # integrate
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)
    # make NES
    nES = iES.mul(intObj.icpVec(), 0)
    nES = np.transpose(nES)
    nES.index = gesObj.obs.index
    
    return(nES)

def Pipeline_Pyther(gesObj, intObj):
    # aREA takes gesObj.X
    nesMat = aREA(gesObj, intObj)

    # Create obs dataframe
    nesMat_sampleNames = pd.DataFrame(index=range(len(nesMat.index.values)),columns=range(0))
    nesMat_sampleNames.index = nesMat.index.values
    nesMat_sampleNames

    # Create var dataframe
    nesMat_proteins = pd.DataFrame(index=range(len(nesMat.columns.values)),columns=range(0))
    nesMat_proteins.index = nesMat.columns.values
    nesMat_proteins

    # Convert the pandas dataframe from Pyther into a new Anndata object
    pAct_adata = anndata.AnnData(X=nesMat,
                                 obs=nesMat_sampleNames,
                                 var=nesMat_proteins)

    # Store the GExpr Anndata object in the PAct Anndata object
    pAct_adata.gExpr = gesObj
    return(pAct_adata)
