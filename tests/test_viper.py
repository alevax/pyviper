
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
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            enrichment="area",
            eset_filter=True
        )
        expected_activity = ad.read_csv(join(rootdir, "test/unit_test_1/viper_area_nes_R.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df())

    def test_viper_narnea(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            enrichment="narnea",
            eset_filter=False
        )
        activity_expected_nes = ad.read_csv(join(rootdir, "test/unit_test_1/viper_narnea_nes_R.csv")).T
        compare_dataframes(activity.to_df(), activity_expected_nes.to_df())

        activity_expected_pes = ad.read_csv(join(rootdir, "test/unit_test_1/viper_narnea_pes_R.csv")).T
        compare_dataframes(activity.to_df(layer="pes"), activity_expected_pes.to_df())

def compare_dataframes(df1, df2, tol=1e-2):
    # Check if the index sets are the same
    if set(df1.index) != set(df2.index):
        raise ValueError("The index sets of the two DataFrames are not the same.")
    
    # Check if the column sets are the same
    if set(df1.columns) != set(df2.columns):
        raise ValueError("The column sets of the two DataFrames are not the same.")
    
    # Align df2 to the order of df1
    df2 = df2.loc[df1.index, df1.columns]
    
    # Compute discrepancies
    differences = np.abs(df1.values - df2.values)
    mean_discrepancy = differences.mean()
    max_discrepancy = differences.max()
    
    # Report discrepancies
    print(f"Mean discrepancy: {mean_discrepancy}")
    print(f"Maximum absolute (MA) discrepancy: {max_discrepancy}")

    # Check if the values are close within tolerance
    if not np.allclose(df1.values, df2.values, atol=tol):
        raise ValueError(f"Numerical values in the DataFrames are not equal within tolerance {tol}.")
    
    print("The two DataFrames are equal within the given tolerance.")

if __name__ == "__main__":
    unittest.main()
