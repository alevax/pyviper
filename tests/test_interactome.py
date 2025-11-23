
import numpy as np
from os.path import dirname, abspath, join
import unittest
import anndata as ad
import scanpy as sc
import pandas as pd

import pyviper

resources_dir = join(dirname(abspath(__file__)), "resources")

class TestInteractome(unittest.TestCase):
    def setUp(self):
        self.table = pd.read_table(join(resources_dir, "test_net1.tsv"), sep="\t")
        self.network = pyviper.Interactome(name="net", net_table=self.table)
    
    def test_str(self):
        print(str(self.network))

    def test_save(self):
        for path in ("network.tsv", "network.csv", "network.pkl"):
            self.network.save(file_path=path)
            saved_network = pyviper.Interactome("saved_net", net_table=path)
            self.assertTrue(np.all(self.network.get_reg_names() == saved_network.get_reg_names()))
            self.assertTrue(np.all(self.network.get_target_names() == saved_network.get_target_names()))
            self.assertTrue(np.all(self.network.net_table["mor"] == saved_network.net_table["mor"]))
            self.assertTrue(np.all(self.network.net_table["likelihood"] == saved_network.net_table["likelihood"]))

    def test_prune(self):
        self.network.prune(max_targets=50, min_targets=5)
        self.assertTrue(self.network.targets_per_regulon().min() >= 5)
        self.assertTrue(self.network.targets_per_regulon().max() <= 50)

    def test_integrate(self):
        network1 = pyviper.Interactome(
            name="net1", 
            net_table=pd.DataFrame({
                'regulator': ['reg1', 'reg2'], 
                'target': ['t1', 't2'], 
                'mor': [0.5, 0.6], 
                'likelihood': [0.7, 0.8]
            })
        )
        network2 = pyviper.Interactome(
            name="net2",
            net_table=pd.DataFrame({
                'regulator': ['reg1', 'reg3'], 
                'target': ['t1', 't3'], 
                'mor': [0.5, 0.6], 
                'likelihood': [0.7, 0.8]
            })
        )
        network1.integrate(network2, normalize_likelihoods=True)
        self.assertEqual(3, network1.size())

    def test_filter(self):
        regulators = self.network.get_reg_names()
        targets = self.network.get_target_names()
        original_network_size = self.network.size()
        original_mean_regulon_size = self.network.targets_per_regulon().mean()

        self.network.filter_targets(targets_remove=targets[2*len(targets)//3:], verbose=True)
        self.assertTrue(self.network.targets_per_regulon().mean() < original_mean_regulon_size)

        self.network.filter_regulators(regulators_remove=regulators[:len(regulators)//2], verbose=True)
        self.assertTrue(self.network.size() < original_network_size)

    def test_translate(self):
        for format in ("human_ensembl", "mouse_ensembl", "mouse_entrez", "human_symbol", "human_ensembl", "human_entrez"):
            original_network_size = self.network.size()
            original_mean_regulon_size = self.network.targets_per_regulon().mean()

            self.network.translate_regulators(desired_format=format, verbose=True)
            self.network.translate_targets(desired_format=format, verbose=True)

            if format == "human_symbol":
                self.assertTrue(self.network.targets_per_regulon().mean() == original_mean_regulon_size)
                self.assertTrue(self.network.size() == original_network_size)
            else:
                self.assertTrue(self.network.targets_per_regulon().mean() <= original_mean_regulon_size)
                self.assertTrue(self.network.size() <= original_network_size)



if __name__ == "__main__":
    unittest.main()