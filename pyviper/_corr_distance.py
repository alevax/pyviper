import numpy as np
import pandas as pd

# SA_GS_subfunctions.R
def _corr_distance_matrix(data):
    # Equivalent to the following in R: d = sqrt(1 - stats::corr(X))
    # Computing the correlation
    corr_matrix = np.corrcoef(data)
    # Calculating sqrt_one_minus_corr_matrix
    sqrt_one_minus_corr_matrix = np.sqrt(1 - corr_matrix)
    # Ensuring diagonal contains 0 values
    np.fill_diagonal(sqrt_one_minus_corr_matrix, 0)
    return(sqrt_one_minus_corr_matrix)

def __add_row_column_names_to_dist_mat(dist_mat, adata):
    # Turn into a dataframe with row and column names
    df_dist = pd.DataFrame(
        dist_mat,
        columns = adata.obs_names,
        index = adata.obs_names
    )
    return(df_dist)

def _corr_distance(adata,
                   use_reduction=True,
                   reduction_slot="X_pca",
                   key_added="corr_dist",
                   copy = False):
    if copy: adata = adata.copy()

    if isinstance(adata, np.ndarray) or isinstance(adata, pd.DataFrame):
        return _corr_distance_matrix(adata)

    if use_reduction == False:
        # use original features
        d = _corr_distance_matrix(adata.X)
    elif use_reduction == True:
        # use principal components
        X = adata.obsm[reduction_slot]
        d = _corr_distance_matrix(X)
    else:
        raise ValueError("reduction must be logical.")
    d = __add_row_column_names_to_dist_mat(d, adata)
    adata.obsp[key_added] = d

    if copy: return adata

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ---------------------------- ** DISTANCE FUNCS ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def corr_distance(adata,
                  use_reduction=True,
                  reduction_slot="X_pca",
                  key_added="corr_dist",
                  copy = False):
    """\
    A tool for computing a distance matrix based on pearson correlation.

    Parameters
    ----------
    adata
        An anndata object containing a signature in adata.X
    use_reduction (default: True)
        Whether to use a reduction (True) (highly recommended - accurate & much faster)
        or to use the direct matrix (False) for computing distance.
    reduction_slot (default: "X_pca")
        If reduction is TRUE, then specify which slot for the reduction to use.
    key_added (default: "corr_dist")
        Slot in obsp to store the resulting distance matrix.

    Returns
    -------
    Adds fields to the input adata, such that it contains a distance matrix
    stored in adata.obsp[key_added].
    """
    return _corr_distance(adata, use_reduction, reduction_slot, key_added, copy)
