
import os
import warnings
warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
from os.path import dirname, abspath, join
import unittest
import anndata as ad
import scanpy as sc

import pyviper
from test_viper import compare_dataframes

resources_dir = join(dirname(abspath(__file__)), "resources")

class TestPyViper(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_csv(join(resources_dir, "ges.csv")).T
        table = pd.read_table(join(resources_dir, "test_net1.tsv"), sep="\t")
        self.network = pyviper.Interactome(name="net1", net_table=table)
        self.network.filter_targets(self.data.var_names)
        self.activity = pyviper.viper(
            gex_data=self.data,
            interactome=self.network,
            enrichment="area",
            min_targets=0,
            eset_filter=True,
            njobs=1,
            verbose=False
        )
        self.expected_activity = pd.read_csv(join(resources_dir, "viper_nes_R_output_RVIPER.csv"), index_col=0).T
        compare_dataframes(self.activity.to_df(), self.expected_activity, max_tol=None, mean_tol=0.01, prefix="aREA NES, ttest")
    
    def test_viper_similarity_two_sided(self):
        similarity = pyviper.pp.viper_similarity(
            nes=self.activity, 
            nn=None,
            ws=(4.0, 2.0), 
            method="two.sided",
            random_state=0, 
            store_in_adata=True, 
            key_added="viper_similarity"
        )
        expected_similarity = pd.read_csv(join(resources_dir, "viper_similarity_R_output_two_sided.csv"), index_col=0).T
        compare_dataframes(similarity, expected_similarity, max_tol=None, mean_tol=0.01, prefix="VIPER Similarity, two-sided")

    def test_viper_similarity_greater_than(self):
        similarity = pyviper.pp.viper_similarity(
            nes=self.activity, 
            nn=None,
            ws=(4.0, 2.0), 
            method="greater",
            random_state=0, 
            store_in_adata=True, 
            key_added="viper_similarity"
        )
        expected_similarity = pd.read_csv(join(resources_dir, "viper_similarity_R_output_greater.csv"), index_col=0).T
        compare_dataframes(similarity, expected_similarity, max_tol=None, mean_tol=0.01, prefix="VIPER Similarity, greater")

    def test_viper_similarity_less_than(self):
        similarity = pyviper.pp.viper_similarity(
            nes=self.activity, 
            nn=None,
            ws=(4.0, 2.0), 
            method="less",
            random_state=0, 
            store_in_adata=True, 
            key_added="viper_similarity"
        )
        expected_similarity = pd.read_csv(join(resources_dir, "viper_similarity_R_output_less.csv"), index_col=0).T
        compare_dataframes(similarity, expected_similarity, max_tol=None, mean_tol=0.01, prefix="VIPER Similarity, two-sided")

    def test_viper_similarity_50_regulators(self):
        similarity = pyviper.pp.viper_similarity(
            nes=self.activity, 
            nn=50,
            ws=(4.0, 2.0), 
            method="less",
            random_state=0, 
            store_in_adata=True, 
            key_added="viper_similarity"
        )
        expected_similarity = pd.read_csv(join(resources_dir, "viper_similarity_R_output_50.csv"), index_col=0).T
        compare_dataframes(similarity, expected_similarity, max_tol=None, mean_tol=0.01, prefix="VIPER Similarity, nn=50, method=less")


if __name__ == "__main__":
    unittest.main()
