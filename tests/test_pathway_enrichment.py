
import os
import numpy as np
import pandas as pd
from os.path import dirname, abspath, join
import unittest
import anndata as ad
import scanpy as sc

import pyviper

rootdir = dirname(dirname(abspath(__file__)))

class TestPathwayEnrichment(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_csv(join(rootdir, "test/unit_test_1/ges.csv")).T
        table = pd.read_table(join(rootdir, "test/unit_test_1/test_net1.tsv"), sep="\t")
        self.network = pyviper.Interactome(name="network", net_table=table)
        self.network.filter_targets(self.data.var_names)
        
    def test_path_enr(self):
        enrichment = pyviper.tl.path_enr(
            self.data, 
            pathway_interactome=self.network, 
            enrichment="narnea", 
            verbose="False", 
            store_input_data=False
        )
        pvalues = pyviper._pp._nes_to_pval_df(
            enrichment.to_df().mean(axis=0), 
            lower_tail=False, 
            adjust=True
        )
        self.assertTrue((pvalues.between(0,1).all()))
    

if __name__ == "__main__":
    unittest.main()
