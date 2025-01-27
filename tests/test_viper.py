
import os
import numpy as np
import pandas as pd
from os.path import dirname, abspath, join
import unittest
import anndata as ad

import pyviper

rootdir = dirname(dirname(abspath(__file__)))

class TestPyViper(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_csv(join(rootdir, "test/unit_test_1/ges.csv")).T
        table1 = pd.read_table(join(rootdir, "test/unit_test_1/test_net1.tsv"), sep="\t")
        table2 = pd.read_table(join(rootdir, "test/unit_test_1/test_net2.tsv"), sep="\t")
        self.network1 = pyviper.Interactome(name="net1", net_table=table1)
        self.network2 = pyviper.Interactome(name="net2", net_table=table2)
        self.network1.filter_targets(self.data.var_names)
        self.network2.filter_targets(self.data.var_names)
    
    def test_viper_area(self):
        activity = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            enrichment="area",
            eset_filter=True
        )
        expected_activity = ad.read_csv(join(rootdir, "test/unit_test_1/viper_area_nes_R.csv")).T
        compare_anndata(activity, expected_activity)

    def test_viper_narnea(self):
        activity = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            enrichment="narnea",
            eset_filter=False
        )
        expected_activity = ad.read_csv(join(rootdir, "test/unit_test_1/viper_narnea_nes_R.csv")).T
        compare_anndata(activity, expected_activity)


def compare_anndata(adata1, adata2, tol=1e-2):
    if set(adata1.obs_names) != set(adata2.obs_names):
        raise ValueError("The sets of obs_names are not the same.")
    
    if set(adata1.var_names) != set(adata2.var_names):
        raise ValueError("The sets of var_names are not the same.")
    
    adata2 = adata2[adata1.obs_names, adata1.var_names]
    
    # Compute discrepancies
    differences = np.abs(adata1.X - adata2.X)
    mean_discrepancy = differences.mean()
    max_discrepancy = differences.max()
    
    # Report discrepancies
    print(f"Mean discrepancy: {mean_discrepancy}")
    print(f"Maximum absolute (MA) discrepancy: {max_discrepancy}")

    if not np.allclose(adata1.X, adata2.X, atol=tol):
        raise ValueError(f"Numerical values in X are not equal within tolerance {tol}.")
    
    print("The two AnnData objects are equal within the given tolerance.")    


if __name__ == "__main__":
    unittest.main()
