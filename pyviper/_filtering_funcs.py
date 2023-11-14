import pandas as pd
import anndata
from .config import config
from ._load._load_regulators import load_TFs, load_coTFs, load_sig, load_surf

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ----------------------------- FILTERING FUNCTIONS ---------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------

# Python program to illustrate the __intersection
# of two lists in most simple way
def __intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

# def __get_mat_from_anndata(adata, layer, features_indices):
#     if layer is None:
#         adata_array = adata.X
#     else:
#         adata_array = adata.layers[layer]
#
#     mat = pd.DataFrame(adata_array[:,features_indices],
#                        index = list(adata.obs_names),
#                        columns = adata.var_names[features_indices,])
#     return(mat)

def __get_features_indices(adata, features_list):
    true_false_list = pd.Series(adata.var_names).isin(features_list).tolist()
    features_indices = [i for i, x in enumerate(true_false_list) if x]
    return(features_indices)

def __get_features_list(adata, feature_groups=None):
    if(type(feature_groups) is str):
        feature_groups = [feature_groups]

    if feature_groups is None:
        features_list = adata.var_names
    else:
        feature_groups = [x.lower() for x in feature_groups]
        features_list = list()
        if "tfs" in feature_groups or "tf" in feature_groups:
            tfs = load_TFs()
            features_list.extend(tfs)
        if "cotfs" in feature_groups or "cotf" in feature_groups:
            cotfs = load_coTFs()
            features_list.extend(cotfs)
        if "sigs" in feature_groups or "sig" in feature_groups:
            sig = load_sig()
            features_list.extend(sig)
        if "surfs" in feature_groups or "surf" in feature_groups:
            surf = load_surf()
            features_list.extend(surf)
        features_list = __intersection(features_list, list(adata.var_names))

    return(features_list)

def __mat_to_anndata(mat):
    # Helper function for *pyviper* and *path_enr*
    # Create obs dataframe
    mat_sampleNames = pd.DataFrame(index=range(len(mat.index.values)),columns=range(0))
    mat_sampleNames.index = mat.index.values
    mat_sampleNames

    # Create var dataframe
    mat_features = pd.DataFrame(index=range(len(mat.columns.values)),columns=range(0))
    mat_features.index = mat.columns.values
    mat_features

    # Convert the pandas dataframe from Pyviper into a new Anndata object
    pax_data = anndata.AnnData(X=mat,
                               obs=mat_sampleNames,
                               var=mat_features)
    return(pax_data)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

# ------------------------- GET ANNDATA OBJECT FILTERED -----------------------
def _get_anndata_filtered_by_feature_group(adata, layer = None, feature_groups=None): #["TFs", "CoTFs", "sig", "surf"],
    features_list = __get_features_list(adata, feature_groups)
    adata_filt = __get_anndata_filtered_by_feature_list(adata, layer, features_list)
    return(adata_filt)

def __get_anndata_filtered_by_feature_list(adata, layer, features_list):
    adata = adata.copy()
    features_indices = __get_features_indices(adata, features_list)
    if layer is not None:
        adata.X = adata.layers[layer]
    adata = adata[:, features_indices]
    return adata
    # mat = __get_mat_from_anndata(adata, layer, features_indices)
    # print(mat)
    # adata_with_features_only = __mat_to_anndata(mat)
    # return(adata_with_features_only)
