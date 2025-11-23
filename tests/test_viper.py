
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

resources_dir = join(dirname(abspath(__file__)), "resources")

class TestPyViper(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_csv(join(resources_dir, "ges.csv")).T
        table1 = pd.read_table(join(resources_dir, "test_net1.tsv"), sep="\t")
        table2 = pd.read_table(join(resources_dir, "test_net2.tsv"), sep="\t")
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
        expected_activity = ad.read_csv(join(resources_dir, "viper_area_nes_R.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df())

    def test_viper_narnea(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            enrichment="narnea",
            eset_filter=False
        )
        activity_expected_nes = ad.read_csv(join(resources_dir, "viper_narnea_nes_R.csv")).T
        compare_dataframes(activity.to_df(), activity_expected_nes.to_df(), prefix="NaRnEA NES")

        activity_expected_pes = ad.read_csv(join(resources_dir, "viper_narnea_pes_R.csv")).T
        compare_dataframes(activity.to_df(layer="pes"), activity_expected_pes.to_df(), prefix="NaRnEA PES")

    def test_viper_area_scale(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            method="scale",
            enrichment="area",
            eset_filter=True
        )
        expected_activity = ad.read_csv(join(resources_dir, "viper_area_nes_scale_R.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df(), max_tol=None, mean_tol=0.01, prefix="aREA NES, scale")

    def test_viper_area_rank(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            method="rank",
            enrichment="area",
            eset_filter=True
        )
        expected_activity = ad.read_csv(join(resources_dir, "viper_area_nes_rank_R.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df(), prefix="aREA NES, rank")

    def test_viper_area_mad(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            method="mad",
            enrichment="area",
            eset_filter=True
        )
        expected_activity = ad.read_csv(join(resources_dir, "viper_area_nes_mad_R.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df(), max_tol=None, mean_tol=0.01, prefix="aREA NES, mad")

    def test_viper_area_ttest(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            method="ttest",
            enrichment="area",
            eset_filter=True
        )
        expected_activity = ad.read_csv(join(resources_dir, "viper_area_nes_ttest_R.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df(), max_tol=None, mean_tol=0.01, prefix="aREA NES, ttest")

    def test_pleiotropy_correction(self):
        activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=self.network1,
            enrichment="area",
            eset_filter=True,
            pleiotropy=True,
            min_targets=0
        )
        expected_activity = ad.read_csv(join(resources_dir, "viper_nes_pleiotropy_corrected_R_output.csv")).T
        compare_dataframes(activity.to_df(), expected_activity.to_df(), max_tol=None, mean_tol=0.01, prefix="aREA Pleiotropy")


def compare_dataframes(df1, df2, max_tol=1e-2, mean_tol=1e-6, prefix=""):
    # Check if the index sets are the same
    if set(df1.index) != set(df2.index):
        raise ValueError("The index sets of the two DataFrames are not the same.")
    
    # Check if the column sets are the same
    if set(df1.columns) != set(df2.columns):
        raise ValueError("The column sets of the two DataFrames are not the same.")
    
    # Align df2 to the order of df1
    df2 = df2.loc[df1.index, df1.columns]
    
    # Compute discrepancies
    abs_discrepancies = np.abs(df1.values - df2.values)
    mean_discrepancy = abs_discrepancies.mean()
    max_discrepancy = abs_discrepancies.max()
     
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_discrepancies = np.where(
            np.abs(df1.values) == 0,
            0,  # or np.nan, depending on whether you want to ignore or flag these
            abs_discrepancies / np.abs(df1.values)
        )
    mean_err = relative_discrepancies.mean()
    max_err = relative_discrepancies.max()

    print(f"Absolute Discrepancy: {prefix}")
    print(pd.Series(abs_discrepancies.flatten().round(4)).quantile(q=[0.0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995, 0.999, 1.0]))
    print(f"Relative Discrepancy: {prefix}")
    print(pd.Series(relative_discrepancies.flatten().round(4)).quantile(q=[0.0, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995, 0.999, 1.0]))
    
    # Report discrepancies
    print(f"Mean absolute discrepancy: {mean_discrepancy}")
    print(f"Maximum absolute (MA) discrepancy: {max_discrepancy}")
    print(f"Mean relative error: {mean_err}")
    print(f"Maximum relative error: {max_err}")

    # Check if the values are all within the maximum tolerance
    if max_tol is not None and max_discrepancy > max_tol:
        raise ValueError(f"Numerical values in the DataFrames are not equal within tolerance {max_tol}.")

    # Check if mean discrepancy is within stated tolerance
    if mean_tol is not None and mean_discrepancy > mean_tol:
        raise ValueError(f"Mean absolute discrepancy between values is not within tolerance {mean_tol}.")

    print("The two DataFrames are equal within the given tolerance.")

if __name__ == "__main__":
    unittest.main()
