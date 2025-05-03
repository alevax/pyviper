import unittest
import numpy as np
import pandas as pd
from os.path import join, dirname, abspath
import anndata as ad
import scanpy as sc
from scipy.stats import rankdata

import pyviper
from pyviper.pp import rank_norm, _median, _mad_from_R

rootdir = dirname(dirname(abspath(__file__)))

class TestRankNorm(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_text(join(rootdir, "test/unit_tests/test_2/test_2_inputs/LCRN1_gExpr_GES.tsv")).T    
        sc.pp.normalize_total(self.data, inplace=True,target_sum=1e6)
        sc.pp.log1p(self.data)
        sc.pp.highly_variable_genes(self.data, flavor="seurat", n_top_genes=2000, inplace=True)
        sc.tl.pca(self.data, svd_solver='arpack', random_state=0)

    def test_rank_norm(self):
        result_median = rank_norm(
            self.data,
            NUM_FUN=_median,
            DEM_FUN =_mad_from_R,
            layer=None,
            key_added=None,
            copy=True
        ).to_df().values
        result_mean = rank_norm(
            self.data,
            NUM_FUN=np.mean,
            DEM_FUN=np.std,
            layer=None,
            key_added=None,
            copy=True
        ).to_df().values

    def test_top_mrs(self):
        # Add an observation column to simulate clustering labels
        # self.adata.obs["cluster"] = ["A", "A", "B", "B", "C"]
        result = pyviper.tl._find_top_mrs(self.data, N=3, both=True, return_as_df=True)


class TestIntegrate(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_csv(join(rootdir, "test/unit_test_1/ges.csv")).T
        table1 = pd.read_table(join(rootdir, "test/unit_test_1/test_net1.tsv"), sep="\t")
        table2 = pd.read_table(join(rootdir, "test/unit_test_1/test_net2.tsv"), sep="\t")
        self.network1 = pyviper.Interactome(name="net1", net_table=table1)
        self.network2 = pyviper.Interactome(name="net2", net_table=table2)
        self.network1.filter_targets(self.data.var_names)
        self.network2.filter_targets(self.data.var_names)

        self.activity: ad.AnnData = pyviper.viper(
            gex_data=self.data, 
            interactome=[self.network1, self.network2],
            enrichment="area",
            eset_filter=True
        )
        pyviper.tl.pca(self.activity, filter_by_feature_groups=["tfs","cotfs"], zero_center=True,  svd_solver='arpack', random_state=0)
        sc.pp.neighbors(self.activity, metric="correlation", n_neighbors=20, n_pcs=50, random_state=0)
        sc.tl.leiden(self.activity, resolution=0.1, n_iterations=-1, random_state=0)
    
    def test_stouffer(self):
        integrated = pyviper.pp.stouffer(self.activity, "leiden", filter_by_feature_groups=["tfs","cotfs"], compute_pvals=True, return_as_df=True)
        pvals = integrated.filter(like="pval", axis=0)
        self.assertEqual(0, integrated.isna().sum().sum())
        self.assertTrue(np.all((0 <= pvals.values) & (pvals.values <= 1)))

    def test_mwu(self):
        integrated = pyviper.pp.mwu(self.activity, "leiden", filter_by_feature_groups=["tfs","cotfs"], compute_pvals=True, return_as_df=True)
        pvals = integrated.filter(like="pval", axis=0)
        stats = integrated.loc[integrated.index.isin(pvals.index), :]
        u_upper_bound = self.activity.obs["leiden"].value_counts().prod()
        self.assertEqual(0, integrated.isna().sum().sum())
        self.assertTrue(np.all((0 <= stats.values) & (stats.values <= u_upper_bound)))
        self.assertTrue(np.all((pvals.values >= 0) & (pvals.values <= 1)))

    def test_spearman(self):
        integrated = pyviper.pp.spearman(self.activity, groupby="leiden", filter_by_feature_groups=["tfs","cotfs"], compute_pvals=True, return_as_df=True)
        self.assertEqual(0, integrated.isna().sum().sum())
        pvals = integrated.filter(like="pval", axis=0)
        stats = integrated.loc[integrated.index.isin(pvals.index), :]
        self.assertTrue(np.all((0 <= pvals.values) & (pvals.values <= 1)))
        self.assertTrue(np.all((-1 <= stats.values) & (stats.values <= 1)))
        
    def test_mean_diffs(self):
        pyviper.pp._mean_diffs(self.activity, "leiden", filter_by_feature_groups=["tfs","cotfs"], return_as_df=False) 
        integrated = pyviper.pp._mean_diffs(self.activity, "leiden", filter_by_feature_groups=["tfs","cotfs"], return_as_df=True)
        self.assertAlmostEqual(np.abs((integrated.values - self.activity.var.T.values)).mean(), 0, places=1)
    
        
if __name__ == '__main__':
    unittest.main()
