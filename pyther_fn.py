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

def aREA(ges, intObj):
    gesMat = ges.X
    
    # rank transform the GES
    rankMat = rankdata(gesMat, axis = 0)
    
    # find intersecting genes
    targetSet = intObj.get_targetSet()
    obsNames = gesObj.obs_names.to_list()
    intersectGenes = [value for value in targetSet if value in obsNames]
    
    # reduce regulon matrices
    icMat = intObj.icMat().loc[intersectGenes]
    morMat = intObj.morMat().loc[intersectGenes]
    print(icMat.shape)
    print(morMat.shape)
    
    # prepare the 1-tailed / 2-tailed matrices
    gesInds = [obsNames.index(i) for i in intersectGenes]
    ges2T = rankMat / (rankMat.shape[0] + 1)
    ges1T = abs(ges2T - 0.5) * 2
    ges1T = ges1T + (1 - np.max(ges1T))/2
    ges2TQ = norm.ppf(ges2T[gesInds, :])
    ges1TQ = norm.ppf(ges1T[gesInds, :])
    print(ges2TQ)
    print(ges1TQ)
    
    # 2-tail enrichment
    dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
    dES = dES.dot(ges2TQ)
    print(dES.shape)
    print(dES)
    
    # 1-tail enrichemnt
    uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
    uES = uES.dot(ges1TQ)
    print(uES.shape)
    print(uES)
    
    # integrate
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)
    print(iES.shape)
    print(iES)
    # make NES
    nES = iES.mul(intObj.icpVec(), 0)
    
    return(nES)






