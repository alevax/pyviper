
import os
from os.path import dirname, abspath, join
import unittest
import anndata

import pyviper

rootdir = dirname(dirname(abspath(__file__)))

class TestPyViper(unittest.TestCase):

    def test_pyviper(self):
        # Load sample data
        ges = anndata.read_text(join(rootdir, "test/unit_tests/test_1/test_1_inputs/LNCaPWT_gExpr_GES.tsv")).T
        pyviper.pp.translate(ges, desired_format="human_symbol")

        # Load network
        network = pyviper.load.msigdb_regulon("h")

        # Filter targets in the interactome
        network.filter_targets(ges.var_names)

        activity = pyviper.viper(gex_data=ges, interactome=network, enrichment="area")

        # narnea
        activity = pyviper.viper(gex_data=ges, interactome=network, enrichment="narnea", eset_filter=False)
        

# Run the tests
if __name__ == "__main__":
    unittest.main()
