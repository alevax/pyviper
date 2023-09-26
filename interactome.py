### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import os

### ---------- EXPORT LIST ----------
__all__ = ['Interactome']

class Interactome:
    # class initialization
    def __init__(self, name, net_table=None, input_type=None):
        self.name = name
        if net_table is None:
            self.net_table = pd.DataFrame(columns=["regulator", "target", "mor", "likelihood"])
        elif type(net_table) is str:
            file_path = net_table
            if input_type is None:
                file_extension = os.path.splitext(file_path)[-1].lower()
                if file_extension == ".csv":
                    self.net_table = pd.read_csv(file_path, sep=",")
                elif file_extension == ".tsv" or file_extension == ".net":
                    self.net_table = pd.read_csv(file_path, sep="\t")
                elif file_extension == ".pkl":
                    self.net_table = pd.read_pickle(file_path)
                else:
                    raise ValueError("Unsupported file format: {}".format(file_extension))
            else:
                if "csv" in input_type:
                    self.net_table = pd.read_csv(file_path, sep=",")
                elif "tsv" in input_type:
                    self.net_table = pd.read_csv(file_path, sep="\t")
                elif "pkl" in input_type or "pickle" in input_type:
                    self.net_table = pd.read_pickle(file_path)
                else:
                    raise ValueError("Unsupported file format: {}".format(input_type))
        else:
            self.net_table = net_table

    # string representation
    def __str__(self):
        retStr = "Object of class Interactome:\n"
        retStr += "\tName: " + self.name
        retStr += "\tNumber of Regulons: " + str(self.size())
        return retStr

    def save(self, file_path, output_type=None):
        if output_type is None:
            file_extension = os.path.splitext(file_path)[-1].lower()
            if file_extension == ".csv":
                self.net_table.to_csv(file_path, sep=",", index=False)
            elif file_extension == ".tsv":
                self.net_table.to_csv(file_path, sep="\t", index=False)
            elif file_extension == ".pkl":
                self.net_table.to_pickle(file_path)
            else:
                raise ValueError("Unsupported file format: {}".format(file_extension))
        else:
            if "csv" in output_type:
                self.net_table.to_csv(file_path, sep=",", index=False)
            elif "tsv" in output_type:
                self.net_table.to_csv(file_path, sep="\t", index=False)
            elif "pkl" in output_type or "pickle" in output_type:
                self.net_table.to_pickle(file_path)
            else:
                raise ValueError("Unsupported file format: {}".format(output_type))

    def copy(self):
        return Interactome(self.name, self.net_table.copy())

    # returns size as the number of regulons
    def size(self):
        return len(self.get_regulonNames())

    def get_regulonNames(self):
        return self.net_table["regulator"].unique()

    def integrate(self, network_list, network_weights = None, normalize_likelihoods = False):
        # self.net_table = pd.concat(self.net_table, net_table)
        if isinstance(network_list, Interactome):
            network_list = [network_list]
        n_networks = len(network_list)
        if network_weights is not None:
            if len(network_weights) != n_networks:
                raise ValueError("network_weights is length" + str(len(network_list)) + ". Should be equal to number of networks: " + str(n_networks) + ".")

        for i in range(0, n_networks):
            if isinstance(network_list[i], Interactome):
                network_list[i] = network_list[i].net_table
            elif not isinstance(network_list[i], pd.DataFrame):
                raise ValueError("Unsupported type of network input:" + str(type(network_list[i])))
        network_list.append(self.net_table)
        n_networks = len(network_list)

        # Get all pairs of regulators and targets
        all_pairs = np.vstack([df[['regulator', 'target']].values for df in network_list])
        unique_pairs_df = pd.DataFrame(all_pairs).drop_duplicates()
        unique_pairs_df.columns = ["regulator", "target"]
        all_pairs_df = unique_pairs_df.sort_values(by=['regulator', 'target'])

        # Make weights if not given
        if network_weights is None:
            network_weights = np.ones(n_networks)

        # Modify each network to have the same regulator-target pairs by
        # adding empty rows of 0s so they be all lined up together in a stack
        net_array_list = []

        for i in range(0, n_networks):
            # Get the missing pairs from a single dataframe
            net_df = network_list[i]
            pairs_in_single_df = net_df[['regulator', 'target']]
            # Convert the columns to sets for faster set operations
            all_pairs_set = set(map(tuple, all_pairs_df[['regulator', 'target']].values))
            single_pairs_set = set(map(tuple, pairs_in_single_df[['regulator', 'target']].values))
            # Find missing pairs using set difference
            missing_pairs_set = all_pairs_set - single_pairs_set
            # Convert the missing pairs back to a DataFrame
            missing_pairs_df = pd.DataFrame(list(missing_pairs_set), columns=['regulator', 'target'])
            # Check if regulators in missing_pairs_df are in net_df["regulator"]
            in_net_df = missing_pairs_df['regulator'].isin(net_df['regulator'])
            # Create two DataFrames based on the condition
            missing_pairs_shared_regs_df = missing_pairs_df[in_net_df]
            missing_pairs_nonshared_regs_df = missing_pairs_df[~in_net_df]

            # If the regulator of the missing pairs is in this network
            # then we want all targets to be set to 0
            missing_pairs_shared_regs_df.insert(2, 'mor', 0)
            missing_pairs_shared_regs_df.insert(3, 'likelihood', 0)

            # If the regulator is not in this network, then this network
            # provide no information and we don't want it affecting means
            # so it to NaN. If it were 0, then it would decrease mor and likelihood
            # of all regulon targets not included with 0s. We don't want that.
            missing_pairs_nonshared_regs_df.insert(2, 'mor', np.nan)
            missing_pairs_nonshared_regs_df.insert(3, 'likelihood', np.nan)

            # Add missing pairs
            net_df = pd.concat([net_df,
                                missing_pairs_shared_regs_df,
                                missing_pairs_nonshared_regs_df], ignore_index=True)
            # Sort so all DataFrames have the same columns for regulator and target
            net_df = net_df.sort_values(by=['regulator', 'target'])
            values_df = net_df.loc[:, ['mor', 'likelihood']].values
            # Multiply by the network weights
            values_df = values_df * network_weights[i]
            net_array_list.append(values_df)

        # Create a stack of networks & compute means across the stack
        network_stack = np.stack(net_array_list)
        stack_means = np.nanmean(network_stack, axis=0)

        # Combine the mean
        df = pd.DataFrame({
            "regulator":all_pairs_df['regulator'],
            "target":all_pairs_df['target'],
            'mor':stack_means[:, 0],
            'likelihood':stack_means[:, 1]
        })

        if normalize_likelihoods:
            # Group by "regulator" and compute ranks of "likelihood"
            df['rank'] = df.groupby('regulator')['likelihood'].rank(ascending=False)

            # Compute the total number of targets + 1 for each regulator
            df['total_targets'] = df.groupby('regulator')['target'].transform('count') + 1

            # Divide "likelihood" by (total number of targets + 1) for each regulator
            df['likelihood'] = df['rank'] / df['total_targets']

            # Drop the intermediate columns if needed
            df = df.drop(['rank', 'total_targets'], axis=1)

        self.net_table = df

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

    def filter_regulators(self, regulators_keep = None, regulators_remove = None):
        if regulators_keep is not None:
            self.net_table = self.net_table[self.net_table['regulator'].isin(regulators_keep)]
        if regulators_remove is not None:
            self.net_table = self.net_table[~self.net_table['regulator'].isin(regulators_remove)]

    def filter_targets(self, targets_keep = None, targets_remove = None):
        if targets_keep is not None:
            self.net_table = self.net_table[self.net_table['target'].isin(targets_keep)]
        if targets_remove is not None:
            self.net_table = self.net_table[~self.net_table['target'].isin(targets_remove)]

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
