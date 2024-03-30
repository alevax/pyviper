import numpy as np
import pandas as pd
from scipy.stats import rankdata
import random
from tqdm import tqdm
from math import log2
from ._corr_distance import corr_distance

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ---------------------------- ** KNN ARRAY FUNC ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-

def __rank_column(column): return rankdata(column, method='ordinal')

def __get_top_n_indices(dist_array_ranked, knn, njobs = 1):
    # Sort each row of the input array and get the indices of the sorted elements
    if njobs == 1:
        sorted_indices = np.argsort(dist_array_ranked, axis=1)
    else:
        sorted_indices = np.array(Parallel(n_jobs=njobs)(
            delayed(np.argsort)(row) for row in dist_array_ranked
        ))

    # Slice the sorted indices to keep only the top knn indices for each row
    top_n_indices = sorted_indices[:, :knn]

    return top_n_indices

def _get_knn_array(dist_df, max_knn, njobs = 1):
    # For each sample, sort neighbors by distance.
    if njobs ==1:
        dist_array_ranked = np.apply_along_axis(
            __rank_column,
            axis=1,
            arr=np.array(dist_df)
        )
    else:
        dist_array_ranked = np.array(Parallel(n_jobs=njobs)(
            delayed(__rank_column)(column) for column in dist_df
        ))

    # Trim this ranking to MaxNN
    max_knn_sorted_indices_array = __get_top_n_indices(dist_array_ranked, max_knn, njobs)

    return max_knn_sorted_indices_array

# ------------------------------ ** HELPER FUNC ** -----------------------------
def _neighbors_knn(adata,
                   max_knn=101,
                   dist_slot="corr_dist",
                   key_added="knn",
                   njobs = 1):
    if isinstance(adata, np.ndarray) or isinstance(adata, pd.DataFrame):
        return _get_knn_array(dist_df = adata, max_knn = max_knn, njobs = njobs)

    # Get the distance DataFrame
    dist_df = adata.obsp[dist_slot]

    # Comptue the KNN array
    max_knn_sorted_indices_array = _get_knn_array(dist_df, max_knn, njobs)

    # Add this to the adata
    adata.uns[key_added] = max_knn_sorted_indices_array

# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
# ------------------------------------------------------------------------------
# ---------------------------- ** KNN ARRAY FUNC ** ----------------------------
# ------------------------------------------------------------------------------
# @-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-@-
def neighbors_knn(adata,
                  max_knn=101,
                  dist_slot="corr_dist",
                  key_added="knn",
                  njobs = 1):
    """\
    A tool for computing a KNN array used to then rapidly generate connectivity
    graphs with acdc.pp.neighbors_graph for clustering.

    Parameters
    ----------
    adata
        An anndata object containing a distance object in adata.obsp.
    max_knn (default: 101)
        The maximum number of k-nearest neighbors (knn) to include in this array.
        acdc.pp.neighbors_graph will only be able to compute KNN graphs with
        knn <= max_knn.
    dist_slot (default: "corr_dist")
        The slot in adata.obsp where the distance object is stored. One way of
        generating this object is with adata.pp.corr_distance.
    key_added (default: "knn")
        Slot in adata.uns to store the resulting knn array.
    njobs (default: 1)
        Paralleization option that allows users to speed up runtime.

    Returns
    -------
    Adds fields to the input adata, such that it contains a knn array stored in
    adata.uns[key_added].
    """
    if isinstance(adata, np.ndarray) or isinstance(adata, pd.DataFrame):
        return _neighbors_knn(adata, max_knn, dist_slot, key_added, njobs)
    else:
        _neighbors_knn(adata, max_knn, dist_slot, key_added, njobs)






def _representative_subsample_pca(pca_array, size, njobs = 1):
    dist_mat = corr_distance(pca_array)
    return representative_subsample_dist(dist_mat, size, njobs)

def representative_subsample_dist(dist_mat, size, njobs = 1):
    n_total_samples = dist_mat.shape[0]
    NN_array = neighbors_knn(dist_mat, n_total_samples, njobs)

    # If sample size is greater than half the number of samples, we take a
    # representative subsample of size (total # cells - sample size) and then
    # remove that representative subsample from all the possible samples. What's
    # left is our final representative subsample that we return.
    if size > int(n_total_samples/2):
        sample_size = n_total_samples - size
    else:
        sample_size = size

    rep_MC_sampling_df = more_equal_sampling_with_nn(NN_array, sample_size, seed = 0)
    first_MC_sample_representative_indices = rep_MC_sampling_df[rep_MC_sampling_df['sample_index'] != -1]['sample_index'].tolist()

    if size > int(n_total_samples/2):
        all_indices = list(range(0, n_total_samples))
        remaining_MC_indices = list(np.setdiff1d(all_indices, first_MC_sample_representative_indices))
        repr_sample = remaining_MC_indices
    else:
        repr_sample = first_MC_sample_representative_indices

    return repr_sample

def more_equal_sampling_with_nn(NN_array, size, seed = 0):
    # We select metacells in a representative way: after selecting each metacell, we eliminate its nearest
    # unselected neighbor from the sampling pool. This way, if you sample more from one subpopulation/cluster of
    # metacells, you are less likely to do so and therefore more likely to sample from the other subpopulations.
    # The reason we use this approach and not just random.sample is because if we have a small population,
    # random.sample may miss it, while using more_equal_sampling_with_nn will capture it.
    n_cells = NN_array.shape[0]
    if size > n_cells:
        raise ValueError("size of " + str(size) + " is larger than the number of samples in NN_array, " + str(n_cells))
    df = pd.DataFrame({'sample_index': [-1] * size, 'nn_index': [-1] * size})
    random.seed(seed)
    randomized_indices = random.sample(list(range(0, n_cells)), n_cells)
    # indices_in_sample = np.array(np.zeros(size), dtype = int)
    excluded_from_sampling = list()
    n_added_samps = 0 # n_collected_samples
    for i in range(0, n_cells):
        index_i = randomized_indices[i]
        if not index_i in excluded_from_sampling:
            # Add index_i as our newest sample
            df.loc[n_added_samps]["sample_index"] = index_i
            # Prevent index_i from being treated as a neighbor
            # of a future sample
            excluded_from_sampling.append(index_i)
            # Identify the NN of index_i that has not
            # yet been added to excluded_from_sampling
            index_i_NN_sorted = NN_array[index_i,:]
            NN_not_excluded = index_i_NN_sorted[~np.isin(index_i_NN_sorted, excluded_from_sampling)]
            first_NN_not_excluded = NN_not_excluded[0]
            df.loc[n_added_samps]["nn_index"] = first_NN_not_excluded
            # Prevent index_i from being treated as a sample or
            # a different sample's neighbor
            excluded_from_sampling.append(first_NN_not_excluded)
            # If we've collected enough samples return our results
            n_added_samps += 1
            if n_added_samps == size:
                break
    return(df)

def condense_in_half(pca_array, metacell_indices_df, seed = 0, njobs = 1):
    pca_array = pca_array.copy()
    n_samples = pca_array.shape[0]
    if n_samples % 2 == 1: # check if odd number of samples
        # Generate a random index for the row to remove
        random.seed(seed)
        random_index = random.randint(0, pca_array.shape[0] - 1)
        # Remove the random row, so we now have an even number of samples
        pca_array = np.delete(pca_array, random_index, axis=0)
        metacell_indices_df = pd.DataFrame(np.delete(metacell_indices_df.values, random_index, axis=0))
        n_samples -= 1

    dist_array = corr_distance(pca_array)

    # Slowest step
    NN_array = neighbors_knn(dist_array, n_samples, njobs)

    n_half_samples = int(n_samples/2)
    sampled_indices_df = more_equal_sampling_with_nn(NN_array, size = n_half_samples, seed = seed)

    # Extract sample_index and nn_index as NumPy arrays
    sample_indices = sampled_indices_df['sample_index'].to_numpy()
    nn_indices = sampled_indices_df['nn_index'].to_numpy()

    # Use NumPy to index rows in pca_array and calculate the mean
    pca_array = np.mean(pca_array[np.vstack((sample_indices, nn_indices))], axis=0)

    # Use NumPy to extract the corresponding rows from metacell_indices_df
    sample_rows = metacell_indices_df.iloc[sample_indices].to_numpy()
    nn_rows = metacell_indices_df.iloc[nn_indices].to_numpy()
    # Stack the rows side by side to create a new array with the
    # original indices in each row, so we can later create the metacells
    metacell_indices_df = pd.DataFrame(np.hstack((sample_rows, nn_rows)))

    return {"pca_array":pca_array, "metacell_indices_df":metacell_indices_df}

def get_mc_indices_by_condensing(pca_array, size = 1000, exact_size = False, seed = 0, verbose = True, njobs = 1): #mode = "exact"
    n_samples = pca_array.shape[0]
    metacell_indices_df = pd.DataFrame(list(range(0, n_samples)))
    # Calculate the number of iterations necessary
    # n_samples/2^n = size ==> n = log2(n_samples / size)
    from math import log2
    total_iters = int(log2(n_samples / size)) # we floor because we want > size
    if verbose: progress_bar = tqdm(total=total_iters, desc="Condensing")
    for _ in range(total_iters): #i.e. while n_samples > size
        condense_results = condense_in_half(pca_array, metacell_indices_df, seed, njobs)
        pca_array = condense_results["pca_array"]
        metacell_indices_df = condense_results["metacell_indices_df"]
        if verbose: progress_bar.update(1)

    if exact_size:
        repr_MC_sample_indices = _representative_subsample_pca(pca_array, size, njobs = 1)
        metacell_indices_df = metacell_indices_df.loc[repr_MC_sample_indices].reset_index(drop=True)

    if verbose: progress_bar.close()

    return metacell_indices_df

# ------------------------------ ** MAIN FUNCS ** ------------------------------

def _representative_subsample_anndata(
    adata,
    pca_slot = "X_pca",
    size = 1000,
    exact_size = True,
    seed = 0,
    key_added = "repr_subsample",
    eliminate = False,
    verbose = True,
    njobs = 1,
    copy = False
):
    if copy: adata = adata.copy()

    pca_array = adata.obsm["X_pca"].copy()
    metacell_indices_df = get_mc_indices_by_condensing(pca_array, size, exact_size, seed, verbose, njobs)
    sample_indices = np.array(metacell_indices_df.loc[:,0])
    adata.obs[key_added] = np.isin(np.arange(adata.shape[0]), sample_indices)
    if eliminate: adata._inplace_subset_obs(sample_indices)
    # adata.uns["knn_groups_indices_df"] = metacell_indices_df

    if copy: return adata
