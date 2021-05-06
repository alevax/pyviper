#!/usr/bin/env python
# coding: utf-8

# In[1]:


from ipynb.fs.full.pyther_classes import *
import numpy as np
from scipy.stats import norm
from scipy.stats import rankdata
import pandas as pd
import anndata 


# In[2]:


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


# In[3]:


def aREA(ges, intObj):
    gesMat = ges.X
    # transform the GES
    ges2t = rankdata(gesMat, axis = 0) / (gesMat.shape[0] + 1)
    ges2tQ = norm.ppf(ges2t)
    ges1t = abs(ges2t - 0.5) * 2
    ges1t = ges1t + (1 - np.max(ges1t)) / 2
    ges1tQ = norm.ppf(ges1t)
    
    # get matrices out of interactome
    icMat = intObj.icMat()
    morMat = intObj.morMat()

    # find the intersecting genes
    # TO DO: flag for not enough genes
    targetSet = intObj.targetSet()
    obsNames = ges.obs_names.to_list()
    intersectGenes = [value for value in targetSet if value in obsNames]
    targetInds = [obsNames.index(i) for i in intersectGenes]
    
    # 2-tail enrichment
    dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
    dES = dES.dot(ges2tQ[targetInds, :])
    # 1-tail enrichemnt
    uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
    uES = uES.dot(ges1tQ[targetInds, :])
    # integrate
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)
    # make NES
    nES = iES.mul(intObj.icpVec(), 0)
    
    return(nES)


# In[4]:


## load data
#intObj = InteractomefromTSV('/mnt/c/Users/lvlah/linux/ac_lab/data/pyther_data/test-net.tsv', 'testNet')
#gesDF = anndata.read_csv('/mnt/c/Users/lvlah/linux/ac_lab/data/pyther_data/gtex-aorta_ges.tsv', delimiter = '\t')
#r = 'ENSG00000000003'


# In[5]:


#aREA(gesDF[:, 0:3], intObj)


# In[ ]:




