import pandas as pd
import numpy as np

class Interactome:
    # class initialization
    def __init__(self, name, netTable=None, sep = '\t'):
        self.name = name
        if netTable is None:
            self.netTable = pd.DataFrame(columns=["regulator", "target", "mor", "likelihood"])
        elif type(netTable) is str:
            self.netTable = pd.read_csv(netTable, sep = sep)
        else:
            self.netTable = netTable

    # string representation
    def __str__(self):
        retStr = "Object of class Interactome:\n"
        retStr += "\tName: " + self.name
        retStr += "\tNumber of Regulons: " + str(self.size())
        return retStr

    def save(self, out_file_path, sep = '\t'):
        self.netTable.to_csv(out_file_path, sep=sep)

    def copy(self):
        return Interactome(self.name, self.netTable.copy())

    # returns size as the number of regulons
    def size(self):
        return len(self.get_regulonNames())

    def get_regulonNames(self):
        return self.netTable["regulator"].unique()

#     def addReg(self, regName, regObj):
#         self.regDict[regName] = regObj

    def get_reg(self, regName):
        return self.netTable[self.netTable['regulator'] == regName]

    # generates the unified set of targets from all regulons
    def get_targetSet(self):
        targetVec = [self.netTable["target"].unique()]
        targetVec = set().union(*targetVec)
        return targetVec


    # generates IC matrix for VIPER
    def icMat(self):
        pivot_df = self.netTable.copy().pivot_table(index='target',
                                                    columns='regulator',
                                                    values='likelihood',
                                                    sort=False,
                                                    fill_value = 0)
        normalization_function = lambda col: col / col.max() / (col / col.max()).sum()
        icMat = pivot_df.apply(normalization_function, axis=0)
        # Reorder columns so they are not alphabetical: keep them consistent with regulonNames: sort=False not working
        icMat = icMat[self.get_regulonNames()]
        return(icMat)

    # generates MoR matrix for VIPER
    def morMat(self):
        morMat = self.netTable.copy().pivot_table(index='target',
                                                  columns='regulator',
                                                  values='mor',
                                                  sort=False,
                                                  fill_value = 0)
        # Reorder columns so they are not alphabetical: keep them consistent with regulonNames: sort=False not working
        morMat = morMat[self.get_regulonNames()]
        return(morMat)

    # generates the vector of icP values for VIPER
    def icpVec(self):
        icP_function = lambda x: np.sqrt(np.sum((x / x.max())**2))
        icP_df = self.netTable.copy().groupby('regulator')['likelihood'].apply(icP_function).reset_index(name='icP')
        unique_regulators = self.get_regulonNames()
        icP_df = icP_df.set_index('regulator').loc[unique_regulators].reset_index()
        icpVec = icP_df["icP"].values
        return icpVec

#     def filter_regulons(self, regulator_list):
#         # First prune the regulators
#         # Using a list comprehension to make a list of the keys to be deleted
#         regulators_to_delete = [key for key in self.regDict if key not in regulator_list]
#         # delete the key/s
#         for regulator in regulators_to_delete:
#             del self.regDict[regulator]

#     def filter_targets(self, targets_list, cleanup = True, min_n_targets = 1):
#         for regulator in self.regDict:
#             regulon = self.regDict[regulator]
#             targets_to_delete = [value for value in regulon.icDict if value not in targets_list]
#             # targets_to_delete = [value for value in ACP1.morDict if value not in targets_list]
#             for target in targets_to_delete:
#                 del regulon.icDict[target]
#                 del regulon.morDict[target]
#         if cleanup is True:
#             # delete the key/s
#             regulators_to_delete = [regulator for regulator in self.regDict if len(self.regDict[regulator].icDict) < min_n_targets]
#             for regulator in regulators_to_delete:
#                 del self.regDict[regulator]

#     # This should be the method from VIPER
#     # later add adaptive and wm parameters
#     def prune(self, cutoff = 50, eliminate = True):
#         for regulator in self.regDict:
#             regulon = self.regDict[regulator]
#             targets_sorted = sorted(regulon.icDict, key=lambda key: regulon.icDict[key], reverse=True)
#             targets_pruned = targets_sorted[:min(len(targets_sorted), cutoff)]
#             regulon.icDict = {target: regulon.icDict[target] for target in targets_pruned}
#             regulon.morDict = {target: regulon.morDict[target] for target in targets_pruned}
#         if eliminate is True:
#             # delete the key/s
#             regulators_to_delete = [regulator for regulator in self.regDict if len(self.regDict[regulator].icDict) < cutoff]
#             for regulator in regulators_to_delete:
#                 del self.regDict[regulator]
