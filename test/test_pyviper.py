
import os
import unittest
import anndata

import pyviper


class TestPyViper(unittest.TestCase):

    def test_pyviper(self):
        # Load sample data
        ges = anndata.read_text("unit_tests/test_1/test_1_inputs/LNCaPWT_gExpr_GES.tsv").T

        # Load network
        network = pyviper.load.msigdb_regulon("h")

        # Translate sample data from ensembl to gene names
        pyviper.pp.translate(ges, desired_format="human_symbol")

        # Filter targets in the interactome
        network.filter_targets(ges.var_names)

        # Compute regulon activities
        # area
        activity = pyviper.viper(gex_data=ges, interactome=network, enrichment="area")

        # narnea
        activity = pyviper.viper(gex_data=ges, interactome=network, enrichment="narnea", eset_filter=False)
        

# Run the tests
if __name__ == "__main__":
    unittest.main()
