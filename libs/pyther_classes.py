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

    def copy(self):
        return Regulon(self.name, self.icDict.copy(), self.morDict.copy())

    # method to create a deep copy of the Regulon object
    def deep_copy(self):
        return Regulon(self.name, copy.deepcopy(self.icDict), copy.deepcopy(self.morDict))
    
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

    def copy(self):
        return Interactome(self.name, self.regDict.copy())

    # method to create a deep copy of the Interactome object
    def deep_copy(self):
        reg_dict_copy = {key: value.deep_copy() for key, value in self.regDict.items()}
        return Interactome(self.name, regDict=reg_dict_copy)
    
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
        targetVec = [r.get_targets() for r in self.regDict.values()]
        targetVec = set().union(*targetVec)
        return(targetVec)

    # generates IC matrix for VIPER
    def icMat(self):
        targetVec = self.get_targetSet()
        # make empty dataframe
        icMat = pd.DataFrame(np.zeros((len(targetVec), self.size())),
                             index = list(targetVec), columns = self.get_regulonNames())
        # fill w/ values from each regulon
        for r in list(self.regDict.values()):
            rNormIC = r.normIC()
            icMat.loc[list(rNormIC.keys()), r.name] = list(rNormIC.values())
        # return matrix
        return(icMat)

    # generates MoR matrix for VIPER
    def morMat(self):
        targetVec = self.get_targetSet()
        # make empty dataframe
        morMat = pd.DataFrame(np.zeros((len(targetVec), self.size())),
                             index = list(targetVec), columns = self.get_regulonNames())
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

    def filter_regulons(self, regulator_list):
        # First prune the regulators
        # Using a list comprehension to make a list of the keys to be deleted
        regulators_to_delete = [key for key in self.regDict if key not in regulator_list]
        # delete the key/s
        for regulator in regulators_to_delete:
            del self.regDict[regulator]

    def filter_targets(self, targets_list, cleanup = True, min_n_targets = 1):
        for regulator in self.regDict:
            regulon = self.regDict[regulator]
            targets_to_delete = [value for value in regulon.icDict if value not in targets_list]
            # targets_to_delete = [value for value in ACP1.morDict if value not in targets_list]
            for target in targets_to_delete:
                del regulon.icDict[target]
                del regulon.morDict[target]
        if cleanup is True:
            # delete the key/s
            regulators_to_delete = [regulator for regulator in self.regDict if len(self.regDict[regulator].icDict) < min_n_targets]
            for regulator in regulators_to_delete:
                del self.regDict[regulator]

    # This should be the method from VIPER
    # later add adaptive and wm parameters
    def prune(self, cutoff = 50, eliminate = True):
        for regulator in self.regDict:
            regulon = self.regDict[regulator]
            targets_sorted = sorted(regulon.icDict, key=lambda key: regulon.icDict[key], reverse=True)
            targets_pruned = targets_sorted[:min(len(targets_sorted), cutoff)]
            regulon.icDict = {target: regulon.icDict[target] for target in targets_pruned}
            regulon.morDict = {target: regulon.morDict[target] for target in targets_pruned}
        if eliminate is True:
            # delete the key/s
            regulators_to_delete = [regulator for regulator in self.regDict if len(self.regDict[regulator].icDict) < cutoff]
            for regulator in regulators_to_delete:
                del self.regDict[regulator]
