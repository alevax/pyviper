#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:34:10 2021

@author: afpvax
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl
# import feather
import igraph
import louvain
import re

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_versions()
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,vector_friendly=False)
# results_file = '/Users/afpvax/Clouds/Dropbox/HDF5/TE001-cells-analysis.h5ad'
data_file = '/Users/afpvax/Clouds/Dropbox/HDF5/TE001-cells-original-data.h5ad'
pd.set_option('display.max_rows', 50)

# Load ADATA
adata = sc.read(data_file)

adata

# ADATA Filtering
sc.pp.filter_cells(adata , min_genes = 500 )
sc.pp.filter_cells(adata , max_genes = 50000 )
sc.pp.filter_cells(adata , min_counts = 1000 )
sc.pp.filter_cells(adata , max_counts = 100000 )
sc.pp.filter_genes(adata , min_cells=5)
adata.shape

# ADATA GEX Analysis
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=9,zero_center=True)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10,metric="euclidean",random_state=666)
sc.tl.louvain(adata,resolution=0.3)
sc.tl.umap(adata,n_components=2,random_state=666)

# Plot UMAP
sc.pl.umap(
    adata,
    #            color=['louvain','Lgr5','Lgr4','Krt19','Olfm4','Dclk1','Bmi1','Gkn3','Hnf4a','Klf5','Cdkn1b','Atoh1'],
    color=[
        "louvain",
        "Lgr5",
        "Lgr4",
        "Krt19",
        "Olfm4",
        "Dclk1",
        "Atoh1",
        "Mki67",
        "Tnfrsf19",
        "Fgfbp1",
    ],
    size=30,
    alpha=0.75,
    ncols=4,
    sort_order=True,
    #            vmax=8 ,
    cmap="viridis",
)


import os
path_str = "/Users/afpvax/Workspace/pyther"
os.chdir(path_str)

from pyther_classes import *
from pyther_fn import *

interactome_filename = "~/Clouds/Dropbox/Data/isc/TE001/networks/TE001-full_unPruned-for-pyther.tsv"

gesObj = adata
# gesObj = anndata.read_csv('ges.tsv', delimiter = '\t')
intObj = InteractomefromTSV( interactome_filename, 'isc_TE001')

## Running VIPER ----
import time
from datetime import timedelta
start_time = time.monotonic()
nesMat = aREA(gesObj, intObj)
end_time = time.monotonic()
print(timedelta(seconds=end_time - start_time))

vp_data = anndata.AnnData(X = nesMat.values ,
                          obs = pd.DataFrame(index=nesMat.index) ,
                          var = pd.DataFrame(index=nesMat.columns) )

sc.tl.pca(vp_data, svd_solver='arpack')
sc.pp.neighbors(vp_data, n_neighbors=15, n_pcs=30,metric="euclidean",random_state=666)
sc.tl.louvain(vp_data,resolution=0.3)
sc.tl.umap(vp_data,n_components=2,random_state=666)

# Plot UMAP
sc.pl.umap(
    vp_data,
    #            color=['louvain','Lgr5','Lgr4','Krt19','Olfm4','Dclk1','Bmi1','Gkn3','Hnf4a','Klf5','Cdkn1b','Atoh1'],
    color=[
        "louvain",
        # "Lgr5",
        # "Lgr4",
        # "Krt19",
        # "Olfm4",
        # "Dclk1",
        "Atoh1",
        # "Mki67",
        # "Tnfrsf19",
        # "Fgfbp1",
    ],
    size=30,
    alpha=0.75,
    ncols=4,
    sort_order=True,
    #            vmax=8 ,
    cmap="viridis",
)


