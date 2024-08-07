# <img width="50" alt="image" src="https://github.com/alevax/pyviper/assets/92543296/3b1f14ab-68a1-49ab-b5b2-aa4f214ad995"> pyVIPER (VIPER Analysis in Python for single-cell RNASeq) 

<!-- [![PyPI](https://img.shields.io/pypi/v/viper-in-python?logo=PyPI)](https://pypi.org/project/viper-in-python) -->
[![PyPI](https://img.shields.io/badge/pypi-stable-darkgreeen?logo=PyPI)](https://pypi.org/project/viper-in-python)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://static.pepy.tech/badge/viper-in-python)](https://pepy.tech/project/viper-in-python)

This package enables network-based protein activity estimation on Python.
It provides also interfaces for scanpy (single-cell RNASeq analysis in Python).
Functions are partly transplanted from R package [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html) and the R package [NaRnEA](https://github.com/califano-lab/NaRnEA). 

The user-friendly documentation is available at: https://alevax.github.io/pyviper/index.html

<div align="center">
    <img style="width: 80%;" alt="viper_visualized" src="https://github.com/user-attachments/assets/44d230cd-7eb3-4cc1-8e27-4cb9d49c269d">
</div>


---

## Dependencies
- `scanpy` for single cell pipeline
- `pandas` and `anndata` for data computing and storage. 
- `numpy` and `scipy`  for scientific computation.
- `joblib` for parallel computing
- `tqdm` show progress bar

If you are using a version of `scanpy` <1.9.3, it is also advisable to downgrade `pandas` to (>=1.3.0 & <2.0), due to `scanpy` incompatibility ([issue](https://github.com/scverse/scanpy/issues/2564))

  
## Installation
### pypi
```shell
pip install viper-in-python
```
### local
```shell
git clone https://github.com/alevax/pyviper/
cd pyviper
pip install -e .
```

## Usage
```python
import pandas as pd
import anndata
import pyviper

# Load sample data
ges = anndata.read_text("test/unit_tests/test_1/test_1_inputs/LNCaPWT_gExpr_GES.tsv").T

# Load network
network = pyviper.load.msigdb_regulon("h")

# Translate sample data from ensembl to gene names
pyviper.pp.translate(ges, desired_format = "human_symbol")

## Filter targets in the interactome
network.filter_targets(ges.var_names)

# Compute regulon activities
## area
activity = pyviper.viper(gex_data=ges, interactome=network, enrichment="area")
print(activity.to_df())

## narnea
activity = pyviper.viper(gex_data=ges, interactome=network, enrichment="narnea", eset_filter=False)
print(activity.to_df())
```

## Tutorials
1. [Analyzing scRNA-seq data at the Protein Activity Level](https://github.com/alevax/pyviper/blob/main/Tutorials/Tutorial-1.ipynb)
2. [Inferring Protein Activity from scRNA-seq data from multiple cell populations with the meta-VIPER approach](https://github.com/alevax/pyviper/blob/main/Tutorials/Tutorial-2.ipynb)
3. [Generating Metacells for ARACNe3 network generation and VIPER protein activity analysis](https://github.com/alevax/pyviper/blob/main/Tutorials/Tutorial-3.ipynb) 

## Structure and rationale

The main functions available from `pyviper` are:
- `pyviper.viper`: "pyviper" function for Virtual Inference of Protein Activity by Enriched Regulon Analysis (VIPER). The function allows using 2 enrichment algorithms, aREA and (matrix)-NaRnEA (see below).
- `pyviper.aREA`: computes [aREA](https://www.nature.com/articles/ng.3593) (analytic rank-based enrichment analysis) and meta-aREA
- `pyviper.NaRnEA`: computes [matrix-NaRnEA](https://www.biorxiv.org/content/10.1101/2021.05.20.445002v5), a vectorized, implementation of [NaRnEA](https://www.mdpi.com/1099-4300/25/3/542)
- `pyviper.pp.translate`: for translating between species (i.e. mouse vs human) and between ensembl, entrez and gene symbols.
- `pyviper.tl.path_enr`: computes pathway enrichment

Other notable functions include:
- `pyviper.tl.OncoMatch`: computes [OncoMatch](https://doi.org/10.1158/2159-8290.CD-22-0342), an algorithm to assess the activity conservation of MR proteins between two sets of samples (e.g. validate GEMMs as effective models of human samples)
- `pyviper.pp.stouffer`: computes signatures on a cluster-by-cluster basis using Cluster integration method for pathway enrichment
- `pyviper.pp.viper_similarity`: computes the [similarity](https://s3.jcloud.sjtu.edu.cn/899a892efef34b1b944a19981040f55b-oss01/bioconductor/3.14/bioc/vignettes/viper/inst/doc/viper.pdf) between VIPER signatures
- `pyviper.pp.repr_metacells`: compute representative metacells (e.g. for ARACNe) using our method to maximize unique sample usage and minimize resampling (users can specify depth, percent data usage, etc).
- `pyviper.pp.repr_subsample`: select a representative subsample of data using our method to ensure a widely distributed sampling.

Additionally, the following submodules are available:
- `pyviper.load`: submodule containing several utility functions useful for different analyses, including `load_msigdb_regulon`, `load_TFs` etc
- `pyviper.pl`: submodule containing pyviper-wrappers for `scanpy` plotting
- `pyviper.tl`: submodule containing pyviper-wrappers for `scanpy` data transformation
- `pyviper.config`: submodule allowing users to specify current species and filepaths for regulators

Last, a new `Interactome` class allows users to load and interrogate ARACNe- and SCENIC-inferred gene regulatory networks.

## Contact
Please, report any issues that you experience through this repository ["Issues"](https://github.com/alevax/pyviper/issues).

For any other info or queries please write to Alessandro Vasciaveo (av2729@cumc.columbia.edu)

## License
`pyviper` is distributed under a MIT License (see [LICENSE](https://github.com/alevax/pyviper/blob/main/LICENSE)).

## Citation
If you used PaCMAP in your publication, or you used the implementation in this repository, please cite our work here:

Zanella, L., Wang, A. L., Lin, Z. J., Vlahos, L. J., Girotto, M. A., Zafar, A., ... & Vasciaveo, A. (2024). PyVIPER: A fast and scalable toolkit for the identification of dysregulated proteins in tumor-derived single-cell RNA-seq data. _Cancer Research, 84_(6_Supplement), 7410-7410.

_Manuscript in review_
