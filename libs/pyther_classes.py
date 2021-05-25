import pandas as pd
import numpy as np

class Regulon:
    # class intialization
    def __init__(self, name, icDict, morDict):
        self.name = name
        self.icDict = icDict
        self.morDict = morDict
    
    # string representation
    def __str__(self):
        retStr = "Object of class Regulon:\n"
        retStr += "\tRegulator: " + self.name
        retStr += "\tNumber of Targets: " + str(len(self.icDict))
        return retStr
    
    # returns size as the number of targets
    def size(self):
        return len(self.icDict.keys())
    
    # returns the set of targets
    def get_targets(self):
        return list(self.icDict.keys())
    
    # returns the normalized wts (VIPER convention) as a dictionary
    def normIC(self):
        # vector or normalized values
        wts = list(self.icDict.values())
        wts = [w / max(wts) for w in wts]
        wts = [w / sum(wts) for w in wts]
        # create dict and return
        wtsDict = dict(zip(self.icDict.keys(), wts))
        return(wtsDict)

    # retuns the likelihood index proportion to all interactions
    def icProportion(self):
        wts = list(self.icDict.values())
        icP = [(w / max(wts))**2 for w in wts]
        icP = sum(icP)**0.5
        return icP
        
class Interactome:
    # class initialization
    def __init__(self, name, regDict=None):
        self.name = name
        if regDict is None:
            self.regDict = {}
        else:
            self.regDict = regDict
    
    # string representation
    def __str__(self):
        retStr = "Object of class Interactome:\n"
        retStr += "\tName: " + self.name
        retStr += "\tNumber of Regulons: " + str(len(self.regDict))
        return retStr
    
    # returns size as the number of regulons
    def size(self):
        return len(self.regDict.keys())

    def get_regulonNames(self):
        return self.regDict.keys()
    
    def addReg(self, regName, regObj):
        self.regDict[regName] = regObj
        
    def get_reg(self, regName):
        return self.regDict[regName]
    
    # generates the unified set of targets from all regulons
    def get_targetSet(self):
        targetVec = [r.targets() for r in self.regDict.values()]
        targetVec = set().union(*targetVec)
        return(targetVec)
    
    # generates IC matrix for VIPER
    def icMat(self):
        targetVec = self.targetSet()
        # make empty dataframe
        icMat = pd.DataFrame(np.zeros((len(targetVec), self.size())), 
                             index = targetVec, columns = self.regulonNames())
        # fill w/ values from each regulon
        for r in list(self.regDict.values()):
            rNormIC = r.normIC()
            icMat.loc[list(rNormIC.keys()), r.name] = list(rNormIC.values())
        # return matrix
        return(icMat)
    
    # generates MoR matrix for VIPER
    def morMat(self):
        targetVec = self.targetSet()
        # make empty dataframe
        morMat = pd.DataFrame(np.zeros((len(targetVec), self.size())), 
                             index = targetVec, columns = self.regulonNames())
        # fill w/ values from each regulon
        for r in list(self.regDict.values()):
            rMOR = r.morDict
            morMat.loc[list(rMOR.keys()), r.name] = list(rMOR.values())
        # return matrix
        return(morMat)
    
    # generates the vector of icP values for VIPER
    def icpVec(self):
        icpVec = [r.icProportion() for r in self.regDict.values()]
        return icpVec


