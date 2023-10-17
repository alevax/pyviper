# Pyther
[![pipy](https://img.shields.io/pypi/v/pyther?color=informational)](https://pypi.python.org/pypi/pyther)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This package enabling network-based protein activity estimation on python. Functions are partly transplanted from R package [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html
).

Find the package user-friendly documentation here: https://califano-lab.github.io/pyther/

---

## Dependencies
- `scanpy` for single cell pipeline
- `pandas` and `anndata` for data computing and storage. 
- `numpy` and `scipy`  for scientific computation.
- `joblib` for parallel computing
- `tqdm` show progress bar

## Installation
### pypi
```shell
pip install pyther
```
### local
```shell
git clone https://github.com/alevax/pyther/tree/main
cd pyther
pip install -e .
```

## Usage
```python
import pandas as pd
import anndata
import pyther

# load sample data
ges = anndata.read_text("test/unit_tests/test_1/test_1_inputs/LNCaPWT_gExpr_GES.tsv").T

# load network
network = pd.read_table("data/regulons/msigdb-h-as-regulon.tsv")
gene_annot = pd.read_csv("data/translate/human_genes.csv", index_col=0)

# prepare interactome
## translate network names from symbol to ensembl
gene_annot = gene_annot[["human_symbol","human_ensembl"]].drop_duplicates().set_index("human_symbol")["human_ensembl"].to_dict()
network["target"] = network["target"].map(gene_annot)
## interactome
network_interactome = pyther.Interactome('net', network)
network_interactome.filter_targets(ges.var_names)

# compute regulon activities
## area
activity = pyther.viper(gex_data=ges, interactome=network_interactome, enrichment="area")
print(activity.to_df())

## narnea
activity = pyther.viper(gex_data=ges, interactome=network_interactome, enrichment="narnea")
print(activity.to_df())
```

## Tutorials
- [Analyzing scRNA-seq data at the Protein Activity Level]()
- [Inferring Protein Activity from scRNA-seq data from multiple cell populations with the meta-VIPER approach]()

## Structure and rationale

The main functions available from `pyther` are:
- `pyther`: "pyther" function for Virtual Inference of Protein Activity by Enriched Regulon Analysis (VIPER). The function allows using 2 enrichment algorithms, aREA and (matrix)-NaRnEA (see below)
- `aREA`: computes [aREA](https://www.nature.com/articles/ng.3593) (analytic rank-based enrichment analysis) and meta-aREA
- `NaRnEA`: computes [matrix-NaRnEA](https://www.biorxiv.org/content/10.1101/2021.05.20.445002v5), a vectorized, implementation of [NaRnEA](https://www.mdpi.com/1099-4300/25/3/542)
- `path_enr`: computes pathway enrichment
- `compute_cluster_stouffer_anndata`: computes signatures on a cluster-by-cluster basis using Cluster integration method for pathway enrichment
- `translate-adata_index`: for mouse-to-human and human-to-mouse

Additionally, the following submodules are available:
- `pyther.load`: submodule containing several utility functions useful for different analyses, including `load_msigdb_regulon`, `load_TFs` etc
- `pyther.pl`: submodule containing pyther-wrappers for `scanpy` plotting
- `pyther.tl`: submodule containing pyther-wrapper for `scanpy` data transformation

## Contact
Please, report any issues that you experience through this repository ["Issues"]().

## License
`pyther` is distributed under a MIT License (see [LICENSE]()).


## Citation


## References
  
