
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
        enrichment = pyviper.tl.path_enrich(
            self.data, 
            interactome=self.network, 
            enrichment="narnea", 
            verbose=False,
            method="ttest"
        )
        pvalues = pyviper._pp._nes_to_pval_df(
            enrichment.to_df().mean(axis=0), 
            lower_tail=False, 
            adjust=True
        )
        self.assertTrue((pvalues.between(0,1).all()))

class TestOncoMatch(unittest.TestCase):

    def test_oncomatch(self):
        # Generate mock data
        rng = np.random.default_rng(42)
        pax_data_to_test = pd.DataFrame(rng.standard_normal((10, 50)), columns=[f'Protein_{i}' for i in range(50)])
        pax_data_for_cMRs = pd.DataFrame(rng.standard_normal((10, 50)), columns=[f'Protein_{i}' for i in range(50)])

        # Run OncoMatch and check output type
        match = pyviper.tl.oncomatch(pax_data_to_test, pax_data_for_cMRs, return_as_df=True)
        self.assertTrue(isinstance(match, pd.DataFrame)),
        self.assertTrue(match.shape == (10, 10))

        # check that each sample is the strongest match for itself
        match = pyviper.tl.oncomatch(pax_data_to_test, pax_data_to_test, return_as_df=True)
        self.assertTrue(np.all(np.diag(match.values) == np.max(match.values, axis=1)))
    

if __name__ == "__main__":
    unittest.main()
