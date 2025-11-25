
import numpy as np
from os.path import dirname, abspath, join
import unittest
import anndata as ad
import scanpy as sc

import pyviper

resources_dir = join(dirname(abspath(__file__)), "resources")

class TestMetacells(unittest.TestCase):
    def setUp(self):
        self.data = ad.read_text(join(resources_dir, "LCRN1_gExpr_GES.tsv")).T
        self.data.raw = self.data
        sc.pp.calculate_qc_metrics(self.data, inplace=True)
        sc.pp.normalize_total(self.data, inplace=True,target_sum=1e6)
        sc.pp.log1p(self.data)
        sc.pp.highly_variable_genes(self.data, flavor="seurat", n_top_genes=2000, inplace=True)
        sc.tl.pca(self.data, svd_solver='arpack', random_state=0)
        pyviper.pp.corr_distance(self.data)
        self.initial_sparsity = np.mean(self.data.raw.X == 0)
        self.initial_mean_depth = np.mean(np.sum(self.data.raw.X, axis=1))

    def test_min_median_depth(self):
        pyviper.pp.repr_metacells(
            adata=self.data,
            counts=None,
            pca_slot="X_pca",
            dist_slot="corr_dist",
            size=10,
            min_median_depth=5000,
            clusters_slot=None,
            key_added=f"metacells",
            seed=12345,
            verbose=False
        )
        metacells = self.data.uns["metacells"]
        metacells_count = metacells.shape[0]
        metacells_sparsity = np.mean(metacells.to_numpy() == 0)
        metacells_n_cells_per = metacells.attrs["n_cells_per_metacell"]

        self.assertTrue(metacells_count, 10)
        self.assertTrue(metacells_sparsity < self.initial_sparsity)
        self.assertEqual(metacells_n_cells_per, int(np.ceil(5000 / self.initial_mean_depth)))

    def test_n_cells_per_metacell(self):
        pyviper.pp.repr_metacells(
            adata=self.data,
            counts=None,
            pca_slot="X_pca",
            dist_slot="corr_dist",
            size=15,
            n_cells_per_metacell=2,
            min_median_depth=None,
            clusters_slot=None,
            key_added=f"metacells",
            seed=12345,
            verbose=False
        )
        metacells = self.data.uns["metacells"]
        metacells_count = metacells.shape[0]
        metacells_sparsity = np.mean(metacells.to_numpy() == 0)
        metacells_median_depth = metacells.attrs["median_depth"]

        self.assertTrue(metacells_count, 15)
        self.assertTrue(metacells_sparsity < self.initial_sparsity)
        self.assertTrue(metacells_median_depth > self.initial_mean_depth * 1.9)

    def test_perc_data_to_use(self):
        pyviper.pp.repr_metacells(
            adata=self.data,
            counts=None,
            pca_slot="X_pca",
            dist_slot="corr_dist",
            size=7,
            n_cells_per_metacell=None,
            perc_data_to_use=100,
            min_median_depth=None,
            clusters_slot=None,
            key_added=f"metacells",
            seed=12345,
            verbose=False
        )
        metacells = self.data.uns["metacells"]
        metacells_count = metacells.shape[0]
        metacells_sparsity = np.mean(metacells.to_numpy() == 0)
        metacells_median_depth = metacells.attrs["median_depth"]

        self.assertTrue(metacells_count, 7)
        self.assertTrue(metacells_sparsity < self.initial_sparsity)
        self.assertTrue(metacells_median_depth > self.initial_mean_depth)

        pyviper.pp.repr_metacells(
            adata=self.data,
            counts=None,
            pca_slot="X_pca",
            dist_slot="corr_dist",
            size=None,
            n_cells_per_metacell=5,
            perc_data_to_use=80,
            min_median_depth=None,
            clusters_slot=None,
            key_added=f"metacells",
            seed=12345,
            verbose=False
        )
        metacells = self.data.uns["metacells"]
        metacells_count = metacells.shape[0]
        metacells_sparsity = np.mean(metacells.to_numpy() == 0)
        metacells_median_depth = metacells.attrs["median_depth"]

        self.assertTrue(metacells_count, 7)
        self.assertTrue(metacells_sparsity < self.initial_sparsity)
        self.assertTrue(metacells_median_depth > self.initial_mean_depth)


    def test_repr_subsample(self):
        sample = pyviper.pp.repr_subsample(
            self.data,
            pca_slot="X_pca",
            size=14,
            seed=0,
            key_added="repr_subsample",
            eliminate=True,
            verbose=True,
            njobs=1,
            copy=True
        )
        self.assertTrue(14, sample.n_obs)
        self.assertTrue(self.data.n_vars, sample.n_vars)


if __name__ == "__main__":
    unittest.main()
