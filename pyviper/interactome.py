### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import os
import pickle
from ._translate import _translate_genes_array
from ._load._load_regulators import _load_regulators

### ---------- EXPORT LIST ----------
__all__ = ['Interactome']

class Interactome:
    # class initialization
    def __init__(self, name, net_table=None, input_type=None):
        """\
        Create an Interactome object to contain the results of ARACNe.
        This object describes the relationship between regulator proteins (e.g.
        TFs and CoTFs) and their downstream target genes with mor (Mode Of
        Regulation, e.g. spearman correlation) indicating directionality and
        likelihood (e.g. mutual information) indicating weight of association.
        An Interactome object can be given to pyviper.viper along with a gene
        expression signature to generate a protein activity matrix with the
        VIPER (Virtual Inference of Protein-activity by Enriched Regulon
        analysis) algorithm[1].

        Parameters
        ----------
        name
            A filepath to one's disk to store the Interactome.
        net_table (default: None)
            Either
            (1) a pd.DataFrame containing four columns in this order:
                "regulator", "target", "mor", "likelihood"
            (2) a filepath to this pd.DataFrame stored either as a .csv,
            .tsv or .pkl.
            (3) a filepath to an Interacome object stored as a .pkl.
        input_type (default: None)
            Only relevant when net_table is a filepath. If None, the input_type
            will be inferred from the net_table. Otherwise, specify "csv", "tsv"
            or "pkl".

        Citations
        -------
        [1] Alvarez, M. J., Shen, Y., Giorgi, F. M., Lachmann, A., Ding, B. B., Ye,
        B. H., & Califano, A. (2016). Functional characterization of somatic
        mutations in cancer using network-based inference of protein activity.
        Nature genetics, 48(8), 838-847.
        """
        self.name = name
        if net_table is None:
            self.net_table = pd.DataFrame(columns=["regulator", "target", "mor", "likelihood"])
        elif type(net_table) is str:
            file_path = net_table
            if input_type is None:
                file_extension = os.path.splitext(file_path)[-1].lower()
                input_type = file_extension[1:]
            if input_type in ["csv", ".csv"]:
                self.net_table = pd.read_csv(file_path, sep=",")
            elif input_type in ["tsv", ".tsv"]:
                self.net_table = pd.read_csv(file_path, sep="\t")
            elif input_type in ["pkl", ".pkl", "pickle"]:
                with open(file_path, 'rb') as file:
                    pkl_obj = pickle.load(file)
                if isinstance(pkl_obj, pd.core.frame.DataFrame):
                    self.net_table = pkl_obj
                elif isinstance(pkl_obj, self. __class__):
                    self.net_table = pkl_obj.net_table
                    self.name = pkl_obj.name
                else:
                    raise ValueError("Unsupported file type: {}".format(type(pkl_obj)))
            else:
                raise ValueError("Unsupported file format: {}".format(input_type))
        else:
            self.net_table = net_table
        n_cols = self.net_table.shape[1]
        if n_cols != 4:
            raise ValueError("net_table contains {} columns".format(n_cols))
        self.net_table.columns = ['regulator', 'target', 'mor', 'likelihood']

    # string representation
    def __str__(self):
        net_table_absMOR = self.net_table.copy()
        net_table_absMOR['mor'] = net_table_absMOR['mor'].abs()
        median_abs_mor_per_reg = str(np.median(net_table_absMOR.groupby('regulator')['mor'].mean()))

        retStr = "Object of class Interactome:"
        retStr += "\n\tName: " + self.name
        retStr += "\n\tNumber of regulons: " + str(self.size())
        retStr += "\n\tMedian average targets per regulon: " + str(np.median(self.net_table.groupby('regulator')['target'].nunique()))
        retStr += "\n\tMedian average abs(mor) per regulon: " + median_abs_mor_per_reg
        retStr += "\n\tMedian average likelihood per regulon: " + str(np.median(self.net_table.groupby('regulator')['likelihood'].mean()))

        return retStr

    def save(self, file_path, output_type=None):
        """\
        Save the Interactome object to one's disk. If saved as "csv" or
        "tsv", just the interactome.net_table will be saved. If saved as
        "pkl", the whole interactome object will be saved.

        Parameters
        ----------
        file_path
            A filepath to one's disk to store the Interactome.
        output_type (default: None)
            If None, the output_type will be inferred from the file_path.
            Otherwise, specify "csv", "tsv" or "pkl".

        Returns
        -------
        None
        """
        if output_type is None:
            file_extension = os.path.splitext(file_path)[-1].lower()
            output_type = file_extension[1:]

        if output_type in ["csv", ".csv"]:
            self.net_table.to_csv(file_path, sep=",", index=False)
        elif output_type in ["tsv", ".tsv"]:
            self.net_table.to_csv(file_path, sep="\t", index=False)
        elif output_type in ["pkl", ".pkl", "pickle"]:
            with open(file_path, 'wb') as file:
                pickle.dump(self, file)
        else:
            raise ValueError("Unsupported file format: {}".format(output_type))

    def copy(self):
        """\
        Create a copy of this Interactome object.

        Returns
        -------
        An object of :class:`~pyviper.interactome.Interactome`.
        """
        return Interactome(self.name, self.net_table.copy())

    # returns size as the number of regulons
    def size(self):
        """\
        Get the the number of regulators in this Interactome.

        Returns
        -------
        An int
        """
        return len(self.get_reg_names())

    def get_reg_names(self):
        """\
        Get an array of all unique regulators in this Interactome.

        Returns
        -------
        An array of strings of :class:`~numpy.ndarray`.
        """
        return self.net_table["regulator"].unique()

    def integrate(self, network_list, network_weights = None, normalize_likelihoods = False):
        """\
        Integrate this Interactome object with one or more other Interacome
        objects to create a consensus network. In general, this should be done
        when interactome objects have the same epigenetics (e.g. due to being
        made from different datasets of same celltype). MetaVIPER should be used
        instead when you have multiple interactomes with different epigenetics
        (e.g. due to being made of data with different celltypes).

        Parameters
        ----------
        network_list
            A single object or a list of objects of class Interactome.
        network_weights (default: None)
            An array containing weights for each network being integrated. The
            first weight corresponds to this network, while the others correspond
            to those in the network list in order. If None, equal weights are
            used.
        normalize_likelihoods (default: False)
            An extra operation that can performed after the integration
            operation where within each regulator, likelihood values are ranked
            and scaled from 0 to 1.
        """
        # self.net_table = pd.concat(self.net_table, net_table)
        if isinstance(network_list, Interactome):
            network_list = [network_list]
        n_networks = len(network_list)

        for i in range(0, n_networks):
            if isinstance(network_list[i], Interactome):
                network_list[i] = network_list[i].net_table
            elif not isinstance(network_list[i], pd.DataFrame):
                raise ValueError("Unsupported type of network input:" + str(type(network_list[i])))
        network_list.insert(0, self.net_table)
        n_networks = len(network_list)

        if network_weights is not None:
            if len(network_weights) != n_networks:
                raise ValueError("network_weights is length" + str(len(network_list)) + ". Should be equal to number of networks: " + str(n_networks) + ".")
            else:
                network_weights = np.array(network_weights)
        else: # Make uniform weights if not given
            network_weights = np.ones(n_networks)

        # Get all pairs of regulators and targets
        all_pairs = np.vstack([df[['regulator', 'target']].values for df in network_list])
        unique_pairs_df = pd.DataFrame(all_pairs).drop_duplicates()
        unique_pairs_df.columns = ["regulator", "target"]
        all_pairs_df = unique_pairs_df.sort_values(by=['regulator', 'target'])

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
            net_array_list.append(values_df)

        # Create a stack of networks & compute means across the stack
        network_stack = np.stack(net_array_list)

        # Recompute the weights as stacks in order to accounts for regulators
        # sets being different in different networks. For each regulator,
        # the weight a network is given will only be counted for when the
        # network contains that regulator.
        weights_stack = network_stack.copy()
        weights_stack[np.isnan(weights_stack) == False] = 1
        weights_stack = weights_stack*network_weights[:, np.newaxis, np.newaxis]
        # For each weight, we normalize to correct to NaN values as follow.
        # We divide by the total across the stacks to get the proportion of
        # each weight to the other, so that they add up to 1. However, we
        # still want the weights to add up to the original total, which they
        # originally didn't with NaN values, so we multiply by the sum of the
        # original weights.
        weights_stack = weights_stack / np.nansum(weights_stack, axis = 0) * np.sum(network_weights)

        # We compute the weighted mean here across the network stacks
        # Multiply each of the items in the stack by their appropriate weight
        network_stack = network_stack * weights_stack
        # Take the sum across the stack and divide by the sum of weights to normalize
        stack_means = np.nansum(network_stack, axis=0)/np.sum(np.array(network_weights))


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
        """\
        Get the rows of the net_table where the regulator is regName.

        Parameters
        ----------
        regName
            The name of a regulator in this Interactome.

        Returns
        -------
        A dataframe of :class:`~pandas.core.frame.DataFrame`.
        """
        return self.net_table[self.net_table['regulator'] == regName]

    # generates the unified set of targets from all regulons
    def get_target_names(self):
        """\
        Get a set of the unique targets in this Interactome

        Returns
        -------
        A 1D NumPy array.
        """
        return self.net_table["target"].unique()

    # generates IC matrix for VIPER
    def ic_mat(self):
        """\
        Get the DataFrame of all the likelihood values. Targets are in the rows,
        while Regulators are in the columns.

        Returns
        -------
        A dataframe of :class:`~pandas.core.frame.DataFrame`.
        """
        pivot_df = self.net_table.copy().pivot_table(index='target',
                                                    columns='regulator',
                                                    values='likelihood',
                                                    sort=False,
                                                    fill_value = 0)
        normalization_function = lambda col: col / col.max() / (col / col.max()).sum()
        ic_mat = pivot_df.apply(normalization_function, axis=0)
        # Reorder columns so they are not alphabetical: keep them consistent with regulonNames: sort=False not working
        ic_mat = ic_mat[self.get_reg_names()]
        return(ic_mat)

    # generates MoR matrix for VIPER
    def mor_mat(self):
        """\
        Get the DataFrame of all the correlation values. Targets are in the rows,
        while Regulators are in the columns.

        Returns
        -------
        A dataframe of :class:`~pandas.core.frame.DataFrame`.
        """
        mor_mat = self.net_table.copy().pivot_table(index='target',
                                                  columns='regulator',
                                                  values='mor',
                                                  sort=False,
                                                  fill_value = 0)
        # Reorder columns so they are not alphabetical: keep them consistent with regulonNames: sort=False not working
        mor_mat = mor_mat[self.get_reg_names()]
        return(mor_mat)

    # generates the vector of icP values for VIPER
    def icp_vec(self):
        """\
        Get the vector containing the proportion of the "Interaction Confidence"
        (IC) score for each interaction in a network, relative to the maximum IC
        score in the network. This vector is generated by taking each individual
        regulon in the newtork and calculating the likelihood index proportion
        to all interactions.

        Returns
        -------
        An array of :class:`~numpy.ndarray`.
        """
        icP_function = lambda x: np.sqrt(np.sum((x / x.max())**2))
        icP_df = self.net_table.copy().groupby('regulator')['likelihood'].apply(icP_function).reset_index(name='icP')
        unique_regulators = self.get_reg_names()
        icP_df = icP_df.set_index('regulator').loc[unique_regulators].reset_index()
        icp_vec = icP_df["icP"].values
        return icp_vec

    def __check_if_reg_names_are_groups(self, regulator_names):
        if len(regulator_names) > 4:
            return False
        else:
            regulator_names = [x.lower() for x in regulator_names]
            all_groups = ["tfs", "cotfs", "sig", "surf"]
            if np.all(np.isin(regulator_names, all_groups)):
                return True
            else:
                return False

    def __get_regulators_by_groups(self, regulator_groups):
        regulator_names = []
        regulator_groups = [x.lower() for x in regulator_groups]
        all_groups = ["tfs", "cotfs", "sig", "surf"]
        if not np.all(np.isin(regulator_names, all_groups)):
            raise ValueError('Given groups are not all one of the following: "tfs", "cotfs", "sig" and "surf".')
        for group in regulator_groups:
            regulator_names = regulator_names + _load_regulators(group)
        return regulator_names


    def filter_regulators(self, regulators_keep = None, regulators_remove = None, verbose = True):
        """\
        Filter regulators by choosing by name or by group which ones you intend
        to keep and which ones you intend to remove from this Interactome.

        Note that the names of regulator that belong to the groups "tfs", "cotfs",
        "sig" and "surf" will be sourced via the paths specified in pyviper.config.
        To update these paths, use the pyviper.config.set_regulators_filepath
        function.

        Parameters
        ----------
        regulators_keep (default: None)
            This should be either:
            (1) An array or list containing the names of specific regulators you
            wish to keep in the network. When left as None, this parameter is
            not used to filter.
            (2) An array or list containing a group or groups of regulators that
            you wish to keep in the network. These groups should be one of the
            following: "tfs", "cotfs", "sig", "surf".
        regulators_remove (default: None)
            This should be either:
            (1) An array or list containing the names of specific regulators you
            wish to remove from the network. When left as None, this parameter
            is not used to filter.
            (2) An array or list containing a group or groups of regulators that
            you wish to remove from the network. These groups should be one of the
            following: "tfs", "cotfs", "sig", "surf".
        verbose (default: True)
            Report the number of regulators removed during filtering
        """
        # Get the number of regulators before filtering
        if verbose:
            n_regs_initial = len(self.net_table["regulator"])
            n_targets_initial = len(self.net_table["target"])

        if regulators_keep is not None:
            if self.__check_if_reg_names_are_groups(regulators_keep) is True:
                regulators_keep = self.__get_regulators_by_groups(regulator_groups = regulators_keep)
            self.net_table = self.net_table[self.net_table['regulator'].isin(regulators_keep)]

        if regulators_remove is not None:
            if self.__check_if_reg_names_are_groups(regulators_remove) is True:
                regulators_remove = self.__get_regulators_by_groups(regulator_groups = regulators_remove)
            self.net_table = self.net_table[~self.net_table['regulator'].isin(regulators_remove)]

        # Report the number of regulators removed during filtering
        if verbose:
            n_targets_final = len(self.net_table["target"])
            n_targets_removed = n_targets_initial - n_targets_final
            print("Removed " + str(n_targets_removed) + " targets.")
            n_regs_final = len(self.net_table["regulator"].unique())
            n_regs_removed = n_regs_initial - n_regs_final
            print("Removed " + str(n_regs_removed) + " regulators.")

    def filter_targets(self, targets_keep = None, targets_remove = None, verbose = True):
        """\
        Filter targets by choosing by name which ones you intend to keep and
        which ones you intend remove from this Interactome.

        When working with an anndata object or a gene expression array, it is
        highly recommended to filter the unPruned network before pruning. This
        is to ensure the pruned network contains a consistent number of targets
        per regulator regulator, all of which exist within gex_data. A regulator
        that has more targets than others will have "boosted" NES scores, such
        that they cannot be compared to those with fewer targets.
        For example, with an anndata object named gex_data, one may is suggested
        to do:
            interactome.filter_targets(gex_data.var_names)

        Parameters
        ----------
        targets_keep (default: None)
            An array containing the names of targets you wish to keep in the
            network. When left as None, this parameter is not used to filter.
        targets_remove (default: None)
            An array containing the names of targets you wish to remove from the
            network. When left as None, this parameter is not used to filter.
        verbose (default: True)
            Report the number of targets removed during filtering
        """

        # Get the number of targets before filtering
        if verbose: n_targets_initial = len(self.net_table["target"])

        if targets_keep is not None:
            self.net_table = self.net_table[self.net_table['target'].isin(targets_keep)]
        if targets_remove is not None:
            self.net_table = self.net_table[~self.net_table['target'].isin(targets_remove)]

        # Report the number of targets removed during filtering
        if verbose:
                n_targets_final = len(self.net_table["target"])
                n_targets_removed = n_targets_initial - n_targets_final
                print("Removed " + str(n_targets_removed) + " targets.")

    def prune(self, max_targets = 50, min_targets = None, eliminate = True, verbose = True):
        """\
        Prune the Interactome by eliminating extra targets from regulators and,
        with eliminate = True, remove regulators with too few targets from the
        network. Note that by ensuring the pruned networks contains the same
        number of targets for each regulator, NES scores are comparable. If one
        regulator has more targest than another, than its NES score will be
        "boosted" and they cannot be compared against each other.

        Parameters
        ----------
        max_targets (default: 50)
            The maximum number of targets that each regulon is allowed.
        min_targets (default: None)
            The minimum number of targets that each regulon is required.
        eliminate (default: True)
            If eliminate = True, then any regulators with fewer targets than
            max_targets will be removed from the network. In other words, after
            pruning, all regulators will have exactly max_targets number of
            targets. This essentially sets min_targets equal to max_targets
            and ensures all NES scores are comparable with aREA.
        verbose (default: True)
            Report the number of targets and regulators removed during pruning
        """
        # Get the number of targets and regulators before pruning
        if verbose:
                n_targets_initial = len(self.net_table["target"])
                n_regs_initial = len(self.net_table["regulator"].unique())

        net_table = self.net_table

        if max_targets is not None:
            # Sort the DataFrame by 'regulator' and 'likelihood' columns
            sorted_df = net_table.sort_values(by=['regulator', 'likelihood'], ascending=[True, False])

            # Group by 'regulator' and apply a function to keep the top 'max_targets' rows in each group
            pruned_df = sorted_df.groupby('regulator').apply(lambda x: x.iloc[:max_targets])

            # Reset the index to flatten the grouped DataFrame
            pruned_df = pruned_df.reset_index(drop=True)

            # Update net_table
            net_table = pruned_df

        if eliminate: min_targets = max_targets
        if min_targets is not None:
            # Count the number of targets for each regulator
            regulator_counts = net_table['regulator'].value_counts()

            # Get the list of regulators with enough targets
            regulators_to_keep = regulator_counts[regulator_counts >= min_targets].index

            # Filter the DataFrame to keep only those regulators
            pruned_df = net_table[net_table['regulator'].isin(regulators_to_keep)]

            # Update net_table
            net_table = pruned_df

        self.net_table = net_table

        # Report the number of targets and regulators after pruning
        if verbose:
                n_targets_final = len(self.net_table["target"])
                n_targets_removed = n_targets_initial - n_targets_final
                print("Removed " + str(n_targets_removed) + " targets.")
                n_regs_final = len(self.net_table["regulator"].unique())
                n_regs_removed = n_regs_initial - n_regs_final
                print("Removed " + str(n_regs_removed) + " regulators.")

    ### ---------- HELPER FUNCTION ----------
    def __translate_net_table_column(self, net_table, desired_format, column_name):
        current_gene_names = net_table[column_name].values
        translation = _translate_genes_array(current_gene_names, desired_format)

        # Update the interactome table with translated values
        net_table[column_name] = translation

        # Removes values without translation
        net_table = net_table.dropna(subset=[column_name])

        return net_table

    def translate_targets(self, desired_format, verbose = True):
        """\
        Translate the targets of the Interactome.  The current name format of
        the targets should be one of the following:
            mouse_symbol, mouse_ensembl, mouse_entrez, human_symbol, human_ensembl or human_entrez

        It is recommended to do this before pruning to ensure a consistent number
        of targets because if targets do not have a translation, they will be
        deleted, resulting in different numbers of targets in a pruned
        interactome that once had consistent number of targets.

        Parameters
        ----------
        desired_format
            Desired format can be one of four strings: "mouse_symbol",
            "mouse_ensembl", "mouse_entrez", "human_symbol", "human_ensembl" or "human_entrez".
        verbose (default: True)
            Report the number of targets successfully and unsucessfully translated
        """
        # Get the number of regulators before translation
        if verbose: n_targets_initial = len(self.net_table["target"].unique())

        self.net_table = self.__translate_net_table_column(
            net_table = self.net_table,
            desired_format = desired_format,
            column_name = 'target'
        )

        # Report how many regulators were removed during translation
        if verbose:
            n_targets_final = len(self.net_table["target"].unique())
            n_targets_removed = n_targets_initial - n_targets_final
            print(str(n_targets_final) + " targets were successfully translated.")
            print(str(n_targets_removed) + " targets were unable to be translated.")

    def translate_regulators(self, desired_format, verbose = True):
        """\
        Translate the regulators of the Interactome. The current name format of
        the regulators should be one of the following:
            mouse_symbol, mouse_ensembl, mouse_entrez, human_symbol, human_ensembl or human_entrez

        Parameters
        ----------
        desired_format
            Desired format can be one of four strings: "mouse_symbol",
            "mouse_ensembl", "mouse_entrez", "human_symbol", "human_ensembl" or "human_entrez".
        verbose (default: True)
            Report the number of regulators successfully and unsucessfully translated
        """
        # Get the number of regulators before translation
        if verbose: n_regs_initial = len(self.net_table["regulator"].unique())

        self.net_table = self.__translate_net_table_column(
            net_table = self.net_table,
            desired_format = desired_format,
            column_name = 'regulator'
        )

        # Report how many regulators were removed during translation
        if verbose:
            n_regs_final = len(self.net_table["regulator"].unique())
            n_regs_removed = n_regs_initial - n_regs_final
            print(str(n_regs_final) + " regulators were successfully translated.")
            print(str(n_regs_removed) + " regulators were unable to be translated.")
