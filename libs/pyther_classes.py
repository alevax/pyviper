import pandas as pd
import numpy as np
import os

class Interactome:
    # class initialization
    def __init__(self, name, net_table=None):
        self.name = name
        if net_table is None:
            self.net_table = pd.DataFrame(columns=["regulator", "target", "mor", "likelihood"])
        elif type(net_table) is str:
            file_path = net_table
            file_extension = os.path.splitext(file_path)[-1].lower()

            if file_extension == ".csv":
                self.net_table = pd.read_csv(file_path, sep=",")
            elif file_extension == ".tsv":
                self.net_table = pd.read_csv(file_path, sep="\t")
            elif file_extension == ".pkl":
                self.net_table = pd.read_pickle(file_path)
            else:
                raise ValueError("Unsupported file format: {}".format(file_extension))
        else:
            self.net_table = net_table

    # string representation
    def __str__(self):
        retStr = "Object of class Interactome:\n"
        retStr += "\tName: " + self.name
        retStr += "\tNumber of Regulons: " + str(self.size())
        return retStr

    def save(self, file_path):
        file_extension = os.path.splitext(file_path)[-1].lower()
        if file_extension == ".csv":
            self.net_table.to_csv(file_path, sep=",", index=False)
        elif file_extension == ".tsv":
            self.net_table.to_csv(file_path, sep="\t", index=False)
        elif file_extension == ".pkl":
            self.net_table.to_pickle(file_path)
        else:
            raise ValueError("Unsupported file format: {}".format(file_extension))

    def copy(self):
        return Interactome(self.name, self.net_table.copy())

    # returns size as the number of regulons
    def size(self):
        return len(self.get_regulonNames())

    def get_regulonNames(self):
        return self.net_table["regulator"].unique()

    def add_regs(self, net_table):
        self.net_table = pd.concat(self.net_table, net_table)

    def get_reg(self, regName):
        return self.net_table[self.net_table['regulator'] == regName]

    # generates the unified set of targets from all regulons
    def get_targetSet(self):
        targetVec = [self.net_table["target"].unique()]
        targetVec = set().union(*targetVec)
        return targetVec


    # generates IC matrix for VIPER
    def icMat(self):
        pivot_df = self.net_table.copy().pivot_table(index='target',
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
        morMat = self.net_table.copy().pivot_table(index='target',
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
        icP_df = self.net_table.copy().groupby('regulator')['likelihood'].apply(icP_function).reset_index(name='icP')
        unique_regulators = self.get_regulonNames()
        icP_df = icP_df.set_index('regulator').loc[unique_regulators].reset_index()
        icpVec = icP_df["icP"].values
        return icpVec

    def filter_regulons(self, regulator_list):
        self.net_table = self.net_table[self.net_table['regulator'].isin(regulator_list)]

    def filter_targets(self, target_list):
        self.net_table = self.net_table[self.net_table['target'].isin(target_list)]

    def prune(self, cutoff = 50, eliminate = True):
        # Sort the DataFrame by 'regulator' and 'likelihood' columns
        sorted_df = self.net_table.sort_values(by=['regulator', 'likelihood'], ascending=[True, False])

        # Group by 'regulator' and apply a function to keep the top 'cutoff' rows in each group
        pruned_df = sorted_df.groupby('regulator').apply(lambda x: x.iloc[:cutoff])

        # Reset the index to flatten the grouped DataFrame
        pruned_df = pruned_df.reset_index(drop=True)

        if eliminate:
            # Count the number of targets for each regulator
            regulator_counts = pruned_df['regulator'].value_counts()

            # Get the list of regulators with enough targets
            regulators_to_keep = regulator_counts[regulator_counts >= cutoff].index

            # Filter the DataFrame to keep only those regulators
            pruned_df = pruned_df[pruned_df['regulator'].isin(regulators_to_keep)]

        self.net_table = pruned_df
