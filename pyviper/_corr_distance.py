import numpy as np
import pandas as pd
from tqdm import tqdm

# SA_GS_subfunctions.R
def _2d_array_vs_2d_array_corr(matrix1, matrix2):
    # Standardize the matrices by subtracting mean and dividing by std (axis=1)
    matrix1_standardized = (matrix1 - matrix1.mean(axis=1, keepdims=True)) / matrix1.std(axis=1, keepdims=True)
    matrix2_standardized = (matrix2 - matrix2.mean(axis=1, keepdims=True)) / matrix2.std(axis=1, keepdims=True)

    # Compute dot product
    dot_product = np.dot(matrix1_standardized, matrix2_standardized.T)

    # Normalize by the number of columns and the vector magnitudes
    magnitudes_matrix1 = np.linalg.norm(matrix1_standardized, axis=1, keepdims=True)
    magnitudes_matrix2 = np.linalg.norm(matrix2_standardized, axis=1, keepdims=True)

    # Compute correlation by dividing the dot product by magnitudes
    correlation_matrix = dot_product / (magnitudes_matrix1 * magnitudes_matrix2.T)

    return correlation_matrix

def _corr_distance_matrix_batch(data, batch_size=1000, verbose=True):
    if isinstance(data, pd.DataFrame):
        data = data.values

    n_samps = data.shape[0]
    sqrt_one_minus_corr_matrix = np.zeros((n_samps, n_samps))

    n_batches = int(np.ceil(n_samps/batch_size))

    for i in tqdm(range(n_batches)) if verbose else range(n_batches):
        # Determine the indices for the current batch
        batch_indices = np.arange(
            i*batch_size,
            min((i+1)*batch_size, n_samps)
        )
        # Get the current batch of data
        batch_data = data[batch_indices, :]

        # Compute correlation of the current batch with all samples
        batch_corr = _2d_array_vs_2d_array_corr(batch_data, data)
        one_min_batch_corr = 1 - batch_corr
        one_min_batch_corr[one_min_batch_corr<0] = 0

        # Compute sqrt(1 - correlation) for the batch
        batch_dist = np.sqrt(one_min_batch_corr)

        # Zero out diagonal elements where the batch is compared with itself
        submatrix = batch_dist[:,batch_indices]
        np.fill_diagonal(submatrix, 0)
        batch_dist[:,batch_indices]=submatrix

        # Store the result in the main distance matrix
        sqrt_one_minus_corr_matrix[batch_indices, :] = batch_dist

        # This ensures correct diagonal filling only where the row and column are from the same batch
        np.fill_diagonal(
            sqrt_one_minus_corr_matrix[np.ix_(batch_indices, batch_indices)], 0
        )

    return sqrt_one_minus_corr_matrix

def _corr_distance_matrix_whole(data):
    # Equivalent to the following in R: d = sqrt(1 - stats::corr(X))
    # Computing the correlation
    corr_matrix = np.corrcoef(data)
    # Calculating sqrt_one_minus_corr_matrix
    sqrt_one_minus_corr_matrix = np.sqrt(1 - corr_matrix)
    # Ensuring diagonal contains 0 values
    np.fill_diagonal(sqrt_one_minus_corr_matrix, 0)
    return(sqrt_one_minus_corr_matrix)

def _corr_distance_matrix(data, batch_size = 1000, verbose = True):
    n_samps = data.shape[0]
    if n_samps < batch_size:
        return _corr_distance_matrix_whole(data)
    else:
        return _corr_distance_matrix_batch(data, batch_size, verbose)

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
                   batch_size=1000,
                   verbose=True,
                   copy = False):
    if copy: adata = adata.copy()
    if isinstance(adata, np.ndarray) or isinstance(adata, pd.DataFrame):
        return _corr_distance_matrix(adata)

    if use_reduction == False:
        # use original features
        d = _corr_distance_matrix(adata.X, batch_size, verbose)
    elif use_reduction == True:
        # use principal components
        X = adata.obsm[reduction_slot]
        d = _corr_distance_matrix(X, batch_size, verbose)
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
                  batch_size=1000,
                  verbose=True,
                  copy=False):
    """\
    A tool for computing a distance matrix based on pearson correlation.

    Parameters
    ----------
    adata
        An anndata object containing a signature in adata.X
    use_reduction : default: True
        Whether to use a reduction (True) (highly recommended - accurate & much faster)
        or to use the direct matrix (False) for computing distance.
    reduction_slot : default: "X_pca"
        If reduction is TRUE, then specify which slot for the reduction to use.
    key_added : default: "corr_dist"
        Slot in obsp to store the resulting distance matrix.
    batch_size : default: 1000
        Reduce total memory usage by running data in batches.
    verbose : default: True
        Show a progress bar for each batch of data.
    copy : default: False
        Return a copy of adata

    Returns
    -------
    Adds fields to the input adata, such that it contains a distance matrix stored in adata.obsp[key_added].
    """
    return _corr_distance(
        adata,
        use_reduction,
        reduction_slot,
        key_added,
        batch_size,
        verbose,
        copy
    )
