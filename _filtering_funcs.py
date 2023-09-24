import pandas as pd

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

def __get_mat_from_anndata(adata, features_indices):
    mat = pd.DataFrame(adata.X[:,features_indices],
                       index = list(adata.obs_names),
                       columns = adata.var_names[features_indices,])
    return(mat)

def __get_features_indices(adata, features_list):
    true_false_list = pd.Series(adata.var_names).isin(features_list).tolist()
    features_indices = [i for i, x in enumerate(true_false_list) if x]
    return(features_indices)

def __get_features_list(adata,
                    feature_groups="all",
                    path_to_tfs = None,
                    path_to_cotfs = None,
                    path_to_sig = None,
                    path_to_surf = None):
    if(type(feature_groups) is str):
        feature_groups = [feature_groups]
    if "all" not in feature_groups:
        feature_groups = [x.lower() for x in feature_groups]
        features_list = list()
        if "tfs" in feature_groups or "tf" in feature_groups:
            tfs = load_TFs(path_to_tfs)
            features_list.extend(tfs)
        if "cotfs" in feature_groups or "cotf" in feature_groups:
            cotfs = load_coTFs(path_to_cotfs)
            features_list.extend(cotfs)
        if "sigs" in feature_groups or "sig" in feature_groups:
            sig = load_sig(path_to_sig)
            features_list.extend(sig)
        if "surfs" in feature_groups or "surf" in feature_groups:
            surf = load_surf(path_to_surf)
            features_list.extend(surf)
        features_list = __intersection(features_list, list(adata.var_names))
    else:
        features_list = adata.var_names
    return(features_list)

def __mat_to_anndata(mat):
    # Helper function for *pyther* and *path_enr*
    # Create obs dataframe
    mat_sampleNames = pd.DataFrame(index=range(len(mat.index.values)),columns=range(0))
    mat_sampleNames.index = mat.index.values
    mat_sampleNames

    # Create var dataframe
    mat_features = pd.DataFrame(index=range(len(mat.columns.values)),columns=range(0))
    mat_features.index = mat.columns.values
    mat_features

    # Convert the pandas dataframe from Pyther into a new Anndata object
    pax_data = anndata.AnnData(X=mat,
                               obs=mat_sampleNames,
                               var=mat_features)
    return(pax_data)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

# ------------------------- GET ANNDATA OBJECT FILTERED -----------------------
def __get_anndata_filtered_by_feature_group(adata,
                               feature_groups="all", #["TFs", "CoTFs", "sig", "surf"],
                               path_to_tfs = None,
                               path_to_cotfs = None,
                               path_to_sig = None,
                               path_to_surf = None):
    features_list = __get_features_list(adata,
            feature_groups,
            path_to_tfs,
            path_to_cotfs,
            path_to_sig,
            path_to_surf)
    adata_filt = __get_anndata_filtered_by_feature_list(adata, features_list)
    return(adata_filt)

def __get_anndata_filtered_by_feature_list(adata, features_list):
    features_indices = __get_features_indices(adata, features_list)
    mat = __get_mat_from_anndata(adata, features_indices)
    adata_with_features_only = __mat_to_anndata(mat)
    return(adata_with_features_only)
