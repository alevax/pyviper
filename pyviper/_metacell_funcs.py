### ---------- IMPORT DEPENDENCIES ----------
import numpy as np
import pandas as pd
import scanpy as sc
from ._rep_subsample_funcs import get_mc_indices_by_condensing, condense_in_half, representative_subsample_dist, neighbors_knn, corr_distance
import warnings
from tqdm import tqdm
import anndata
from scipy.sparse._csr import csr_matrix

def get_counts_as_df(counts, adata):
    if counts is None:
        if adata.raw is None:
            raise ValueError("counts must be given as a parameter or adata.raw must contain counts.")
        else:
            if isinstance(adata.raw.X, np.ndarray):
                counts = pd.DataFrame(adata.raw.X)
            elif isinstance(adata.raw.X, csr_matrix):
                counts = pd.DataFrame(adata.raw.X.toarray())
            else:
                raise ValueError("type(adata.raw.X) is " + str(type(adata.raw.X)) + ". Must be np.ndarray or scipy.sparse._csr.csr_matrix.")
            counts.index = adata.raw.obs_names
            counts.columns = adata.raw.var_names
            counts = counts.loc[adata.obs_names]
    elif isinstance(counts, anndata.AnnData):
        counts = counts.to_df()
    elif isinstance(counts, pd.DataFrame):
        counts = counts
    elif isinstance(counts, np.ndarray):
        counts = pd.DataFrame(counts)
    else:
        raise ValueError("counts must be either anndata.AnnData or pd.DataFrame, or np.ndarray.")


    if adata.shape[0] != counts.shape[0]:
        raise ValueError("Different number of samples between counts (" +  \
                         str(counts.shape[0]) + ") and adata (" + str(adata.shape[0]) + ").")

    return counts

def select_optimal_column_of_knn_groups_indices_df(knn_groups_indices_df, knn_array, progress_bar = True):
    # Each column of the knn groups has a different set of members in the groups
    # E.g. the 0th columns has the 0th member of each group, the 1st column has
    # the 1st member of each group. We check each of these sets and use the one
    # that results in the fewest shared KNN as a quick starting point before
    # optimizing further.
    best_col = 0
    metric_opt = knn_groups_indices_df.size
    for i in tqdm(range(knn_groups_indices_df.shape[1]), desc = "opt selection") if progress_bar else range(knn_groups_indices_df.shape[1]):
        sample_indices = np.array(knn_groups_indices_df.loc[:,i])
        metacell_knn = knn_array[sample_indices,:].copy()
        n_times_sample_used = np.unique(metacell_knn.flatten(), return_counts = True)[1]
        # ave_n_times_sample_used correlates pos with % total samples used and
        # correlates neg with % included samples reused
        ave_n_times_sample_used = np.mean(n_times_sample_used)
        if ave_n_times_sample_used < metric_opt:
            best_col = i
            metric_opt = ave_n_times_sample_used
    sample_indices = np.array(knn_groups_indices_df.loc[:,best_col])
    return sample_indices

def __next_index_with_equal_minimum(arr, index_among_equals = 0):
    min_index = np.argmin(arr)
    min_value = arr[min_index]

    # Find the next occurrence of the minimum value
    indices_of_equal_mins = np.where(arr == min_value)[0]

    if indices_of_equal_mins.size > 1:
        if index_among_equals > indices_of_equal_mins.size - 1:
            index_among_equals = indices_of_equal_mins.size - 1
        min_index = indices_of_equal_mins[index_among_equals]
    return min_index

def optimize_selected_sample_from_each_knn_group(sample_indices, knn_groups_indices_df, knn_array, progress_bar = True):
    n_unique_samples_old = len(np.unique(knn_array[sample_indices,:]))
    randomized_equal_positions = False
    index_among_equal_positions = np.zeros(len(sample_indices), dtype=int)
    for j in tqdm(range(100)) if progress_bar else range(100): #20 iters should be enough to optimize
        for i in range(knn_groups_indices_df.shape[0]):
            # Get all the neighbors of current sample_indices except for group i
            mask = np.arange(len(sample_indices)) != i
            metacell_knn = knn_array[sample_indices[mask],:].copy()

            # Get all the possible neighbor sets within group i
            group_i_indices = knn_groups_indices_df.values[i,:]
            group_i_knn = knn_array[group_i_indices,:].copy()

            # Identify which neighbor set within group i has the fewest shared
            # neighbors with the neighbors of all the other sample_indices
            current_knn_indices_to_use = np.unique(metacell_knn.flatten())
            group_i_index_with_min_shared_knn = \
                np.argmin(np.sum(np.isin(group_i_knn, current_knn_indices_to_use),axis = 1))

            # Switch to the neighbor set with the fewest shared neighbors
            if sample_indices[i] != group_i_indices[group_i_index_with_min_shared_knn]:
                sample_indices[i] = group_i_indices[group_i_index_with_min_shared_knn]
                index_among_equal_positions[i] = 0

        # Break the loop if no changes made
        n_unique_samples = len(np.unique(knn_array[sample_indices,:]))

        n_new_unique_samples = n_unique_samples - n_unique_samples_old
        if n_new_unique_samples == 0:
            if randomized_equal_positions:
                break
            else:
                index_among_equal_positions = index_among_equal_positions + 1
                for i in range(knn_groups_indices_df.shape[0]):
                    # Get all the neighbors of current sample_indices except for group i
                    mask = np.arange(len(sample_indices)) != i
                    metacell_knn = knn_array[sample_indices[mask],:].copy()

                    # Get all the possible neighbor sets within group i
                    group_i_indices = knn_groups_indices_df.values[i,:]
                    group_i_knn = knn_array[group_i_indices,:].copy()

                    # Identify which neighbor set within group i has the fewest shared
                    # neighbors with the neighbors of all the other sample_indices
                    current_knn_indices_to_use = np.unique(metacell_knn.flatten())

                    # MAIN CHANGE HERE FROM ABOVE: __next_index_with_equal_minimum
                    group_i_index_with_min_shared_knn = \
                        __next_index_with_equal_minimum(
                            np.sum(np.isin(group_i_knn, current_knn_indices_to_use),axis = 1),
                            index_among_equal_positions[i]
                        )

                    # Switch to the neighbor set with the fewest shared neighbors
                    sample_indices[i] = group_i_indices[group_i_index_with_min_shared_knn]
                randomized_equal_positions = True
        else:
            randomized_equal_positions = False
            n_unique_samples_old = n_unique_samples
    return sample_indices

def get_sample_indices_with_max_unique_samples(knn_groups_indices_df, knn_array, verbose):
    # Get a good starting point for the optimization
    sample_indices = select_optimal_column_of_knn_groups_indices_df(
        knn_groups_indices_df,
        knn_array,
        progress_bar = verbose
    )
    # Optimize by switching to samples whose shared neighbors don't overlap much
    sample_indices = optimize_selected_sample_from_each_knn_group(
        sample_indices,
        knn_groups_indices_df,
        knn_array,
        progress_bar = verbose
    )
    return sample_indices

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ---------------------------------------------------------------------------
# -------------- PRESET PARAMS METACELL CONSTRUCTION FUNCTIONS --------------
# ---------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_sample_indices_with_preset_params(
    adata,
    counts,
    pca_slot,
    dist_slot,
    size,
    n_cells_per_metacell,
    min_median_depth,
    exact_size = True,
    seed = 0,
    verbose = True,
    njobs = 1,
    smart_sample = True
):
    # if size is None:
    #     size = int(adata.shape[0] / n_cells_per_metacell)
    # elif n_cells_per_metacell is None:
    #     n_cells_per_metacell = int(np.ceil(adata.shape[0] / size))
    if size is None or n_cells_per_metacell is None:
        raise ValueError("get_sample_indices_with_preset_params," + \
                         " but size and n_cells_per_metacell are not both set.")

    pca_array = adata.obsm[pca_slot].copy()
    # Get the KNN to identify the neighbors of each selected index
    if dist_slot not in list(adata.obsp.keys()):
        warnings.warn(dist_slot + " not in adata.obsp. Computing correlation distance...")
        corr_distance(adata,
                      use_reduction=True,
                      reduction_slot=pca_slot,
                      key_added='corr_dist')
        dist_slot = 'corr_dist'

    if smart_sample:
        if verbose: print("Identifying a representative sample from the data...")
        knn_groups_indices_df = get_mc_indices_by_condensing(pca_array, size, exact_size, seed, verbose, njobs)
        if n_cells_per_metacell is None:
            n_cells_per_metacell = knn_groups_indices_df.shape[1]
    else:
        # get random sample using seed
        n_cells_total = adata.shape[0]
        all_indices = np.array(range(n_cells_total))
        np.random.seed(seed)
        sample_indices = np.random.choice(all_indices, size=size, replace=False)
        if n_cells_per_metacell is None:
            n_cells_per_metacell = int(n_cells_total/size)

    if verbose: print("Computing neighbors as knn...")
    neighbors_knn(adata,
                  max_knn=n_cells_per_metacell,
                  dist_slot=dist_slot,
                  key_added='knn4MC',
                  njobs = 1)
    knn_array = adata.uns['knn4MC']

    if verbose: print("Selecting indices for metacells...")
    if smart_sample:
        sample_indices = get_sample_indices_with_max_unique_samples(
            knn_groups_indices_df,
            knn_array,
            verbose
        )
    adata.obs["selected_for_mc"] = pd.Categorical(np.isin(np.arange(adata.shape[0]), sample_indices))
    return adata, size, n_cells_per_metacell

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ---------------------------------------------------------------------------
# ----------- PERCENT OPTIMIZATION METACELL CONSTRUCTION FUNCTIONS ----------
# ---------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ---------------------------------------------------------------------------
# ------- get_sample_indices_by_optimizing_size_by_knn_with_perc_data -------
# ---------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def get_optimally_sized_knn_groups_by_dividing(
    pca_array,
    knn_array,
    perc_data_to_use,
    perc_incl_data_reused,
    seed,
    verbose,
    njobs
):
    n_total_samples = pca_array.shape[0]
    knn_groups_indices_df = pd.DataFrame(list(range(0, n_total_samples)))
    pca_array_original = pca_array.copy()
    sample_indices = None

    max_iters = int(np.ceil(np.log2(n_total_samples)))

    for i in tqdm(range(max_iters)) if verbose else range(max_iters): #prevent while-loop running forever: use large iters for-loop
        condense_results = condense_in_half(pca_array, knn_groups_indices_df, seed, njobs)
        pca_array = condense_results["pca_array"]
        previous_knn_groups_indices_df = knn_groups_indices_df
        knn_groups_indices_df = condense_results["metacell_indices_df"]
        previous_sample_indices = sample_indices
        sample_indices = select_optimal_column_of_knn_groups_indices_df(
            knn_groups_indices_df,
            knn_array,
            progress_bar = False
        )
        sample_indices = optimize_selected_sample_from_each_knn_group(
            sample_indices,
            knn_groups_indices_df,
            knn_array,
            progress_bar = False
        )
        metacell_knn = knn_array[sample_indices,:].copy()
        perc_total_samples_used = len(np.unique(metacell_knn.flatten()))/n_total_samples*100
        # print(i)
        # print(len(sample_indices))
        # print(str(perc_total_samples_used) + "% total samples were used.")

        n_times_sample_used = np.unique(metacell_knn.flatten(), return_counts = True)[1]
        prop_included_samples_reused = np.mean(n_times_sample_used > 1)
        perc_included_samples_reused = np.round(prop_included_samples_reused*100,1)
        # print(str(perc_included_samples_reused) + "% of these included samples were used more than once.")

        if perc_data_to_use is not None:
            if perc_total_samples_used < perc_data_to_use:
                break
        else: #perc_incl_data_reused is not None:
            if perc_included_samples_reused < perc_incl_data_reused:
                break

    return previous_knn_groups_indices_df, previous_sample_indices

def get_sample_indices_by_optimizing_size_by_knn_with_perc_data(
    adata,
    pca_slot = "X_pca",
    dist_slot = "corr_dist",
    n_cells_per_metacell = None,
    perc_data_to_use = None,
    perc_incl_data_reused = None,
    exact_size = True,
    seed = 0,
    verbose = True,
    njobs = 1
):

    pca_array = adata.obsm[pca_slot].copy()

    if dist_slot not in list(adata.obsp.keys()):
        warnings.warn(dist_slot + " not in adata.obsp. Computing correlation distance...")
        corr_distance(adata,
                      use_reduction=True,
                      reduction_slot=pca_slot,
                      key_added='corr_dist')
        dist_slot = 'corr_dist'

    # knn_groups_indices_df = get_mc_indices_by_condensing(pca_array, size, exact_size, seed, verbose, njobs)
    neighbors_knn(adata,
                          max_knn=n_cells_per_metacell,
                          dist_slot=dist_slot,
                          key_added='knn4MC',
                          njobs = 1)
    knn_array = adata.uns['knn4MC']
    previous_sample_indices = None

    # Keep dividing and optimizing the data until you hit the point where you overshot the value
    # e.g. you're aiming for perc_data_to_use = 70. You divide the data and get perc_data_to_use = 80. 80 > 70,
    # so you divide again and get perc_data_to_use = 60. So you stop and take the solution where
    # perc_data_to_use = 80. That solution will then be subsampled (not divided!) until you get
    # perc_data_to_use = 70.
    if verbose: print("Computing initial sample size to optimize...")
    knn_groups_indices_df, upper_bound_potential_samples, = get_optimally_sized_knn_groups_by_dividing(
        pca_array,
        knn_array,
        perc_data_to_use,
        perc_incl_data_reused,
        seed,
        verbose,
        njobs
    )
    n_knn_groups = knn_groups_indices_df.shape[0]
    upper_bound = n_knn_groups # Either use all the samples
    lower_bound = 1 # Or only 1 sample, though that seems unlikely
    sample_indices_potential_dist = adata.obsp['corr_dist'][upper_bound_potential_samples,:][:,upper_bound_potential_samples]

    if verbose: print("Optimizing the sample size to match parameters...")
    for i in tqdm(range(50)) if verbose else range(50):
        est_size = int((upper_bound + lower_bound)/2)

        # print("------------------------")
        # print("est_size")
        # print(est_size)
        # print("upper_bound")
        # print(upper_bound)
        # print("lower_bound")
        # print(lower_bound)

        # Of the KNN groups available that have been optimized by the upper bound,
        # pick a subsample of size est_size
        upper_bound_selected_samples = upper_bound_potential_samples[
            representative_subsample_dist(
                sample_indices_potential_dist,
                est_size
        )]
        # Identify the rows of the KNN groups DF to use based on your subsample
        selected_knn_groups = np.where(np.sum(np.isin(knn_groups_indices_df.values, upper_bound_selected_samples), axis=1) == 1)[0]

        # Identify the optimal sample within each row (KNN group) to use
        temp_sample_indices = select_optimal_column_of_knn_groups_indices_df(
            knn_groups_indices_df.iloc[selected_knn_groups],
            knn_array,
            progress_bar = False
        )
        temp_sample_indices = optimize_selected_sample_from_each_knn_group(
            temp_sample_indices,
            knn_groups_indices_df.iloc[selected_knn_groups],
            knn_array,
            progress_bar = False
        )

        # Create metacells and test their percent sample usage
        metacell_knn = knn_array[temp_sample_indices,:].copy()
        perc_total_samples_used = len(np.unique(metacell_knn.flatten()))/adata.shape[0]*100
        # print(str(perc_total_samples_used) + "% total samples were used.")

        n_times_sample_used = np.unique(metacell_knn.flatten(), return_counts = True)[1]
        prop_included_samples_reused = np.mean(n_times_sample_used > 1)
        perc_included_samples_reused = np.round(prop_included_samples_reused*100,1)
        # print(str(perc_included_samples_reused) + "% of these included samples were used more than once.")

        if perc_data_to_use is not None and\
           round(perc_data_to_use) == int(np.floor(perc_total_samples_used)):
            break
        elif perc_incl_data_reused is not None and\
             round(perc_incl_data_reused) == int(np.ceil(perc_included_samples_reused)):
            break
        elif upper_bound == lower_bound:
            return
            break
        elif int((upper_bound + lower_bound)/2) == lower_bound:
            break
        elif (perc_data_to_use is not None and\
             perc_total_samples_used > perc_data_to_use)\
             or (perc_incl_data_reused is not None and\
             perc_included_samples_reused > perc_incl_data_reused):
            # If we set a new upper bound, narrow down our potential KNN group selection
            upper_bound = est_size
            upper_bound_potential_samples = temp_sample_indices
            sample_indices_potential_dist = adata.obsp['corr_dist'][upper_bound_potential_samples,:][:,upper_bound_potential_samples]
        elif (perc_data_to_use is not None and\
             perc_total_samples_used < perc_data_to_use)\
             or (perc_incl_data_reused is not None and\
             perc_included_samples_reused < perc_incl_data_reused):
            lower_bound = est_size
        else:
            raise ValueError("Error in get_sample_indices_by_optimizing_size_by_knn_with_perc_data")

    adata.obs["selected_for_mc"] = pd.Categorical(np.isin(np.arange(adata.shape[0]), temp_sample_indices))
    size = est_size
    return adata, size, n_cells_per_metacell


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ---------------------------------------------------------------------------
# ------- get_sample_indices_by_optimizing_knn_by_size_with_perc_data -------
# ---------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def optimize_knn_by_size_with_perc_data(
    knn_groups_indices_df,
    knn_array,
    perc_data_to_use,
    perc_incl_data_reused,
    verbose
):
    approx_n_total_samples = knn_groups_indices_df.size
    n_total_samples = knn_array.shape[0]
    if verbose: print("Computing initial knn to optimize...")
    for i in tqdm(range(n_total_samples)) if verbose else range(n_total_samples):
        knn_i_array = knn_array[:,:(i+1)]
        temp_sample_indices = select_optimal_column_of_knn_groups_indices_df(
                knn_groups_indices_df,
                knn_i_array,
                progress_bar = False
            )
        metacell_knn = knn_i_array[temp_sample_indices,:].copy()
        perc_total_samples_used = len(np.unique(metacell_knn.flatten()))/n_total_samples*100
        # print(str(perc_total_samples_used) + "% total samples were used.")

        n_times_sample_used = np.unique(metacell_knn.flatten(), return_counts = True)[1]
        prop_included_samples_reused = np.mean(n_times_sample_used > 1)
        perc_included_samples_reused = np.round(prop_included_samples_reused*100,1)
        # print(str(perc_included_samples_reused) + "% of these included samples were used more than once.")

        if perc_data_to_use is not None and\
           perc_total_samples_used > perc_data_to_use:
            break
        elif perc_incl_data_reused is not None and\
           perc_included_samples_reused > perc_incl_data_reused:
            break

    # SET UP THE BINARY SEARCH ALGORITHM
    if perc_data_to_use is not None:
        # For perc_data_to_use, all values below will be improved b/c of optimization. We will
        # be able to include more data at the same KNN. So we will never go beyond i as the
        # upper bound.
        upper_bound = i
        lower_bound = 1
    elif perc_incl_data_reused is not None:
        # However, for perc_incl_data_reused, it will take many more KNN to reach the
        # same amount of reused data.
        upper_bound = np.min([i * i, n_total_samples])  # mult by 10 b/c above we don't have much optimization
        lower_bound = i

    # Optimize by switching to samples whose shared neighbors don't overlap much

    if verbose: print("Optimizing knn to match parameters...")
    for i in tqdm(range(50)) if verbose else range(50): #max iters to prevent while-loop infinite run
        # BINARY SEARCH:
        est_n_cells_per_metacell = int((upper_bound + lower_bound)/2)
        # print("------------------------")
        # print("est_n_cells_per_metacell")
        # print(est_n_cells_per_metacell)
        # print("upper_bound")
        # print(upper_bound)
        # print("lower_bound")
        # print(lower_bound)
        knn_est_array = knn_array[:,:est_n_cells_per_metacell]
        temp_sample_indices = select_optimal_column_of_knn_groups_indices_df(
                knn_groups_indices_df,
                knn_est_array,
                progress_bar = False
            )
        temp_sample_indices = optimize_selected_sample_from_each_knn_group(
            temp_sample_indices,
            knn_groups_indices_df,
            knn_est_array,
            progress_bar = False
        )
        metacell_knn = knn_est_array[temp_sample_indices,:].copy()
        perc_total_samples_used = len(np.unique(metacell_knn.flatten()))/n_total_samples*100
        # print(str(perc_total_samples_used) + "% total samples were used.")

        n_times_sample_used = np.unique(metacell_knn.flatten(), return_counts = True)[1]
        prop_included_samples_reused = np.mean(n_times_sample_used > 1)
        perc_included_samples_reused = np.round(prop_included_samples_reused*100,1)
        # print(str(perc_included_samples_reused) + "% of these included samples were used more than once.")

        if perc_data_to_use is not None and\
           round(perc_data_to_use) == int(np.floor(perc_total_samples_used)):
            break
        elif perc_incl_data_reused is not None and\
             round(perc_incl_data_reused) == int(np.ceil(perc_included_samples_reused)):
            break
        elif upper_bound == lower_bound:
            break
        elif int((upper_bound + lower_bound)/2) == lower_bound:
            break
        elif (perc_data_to_use is not None and\
             perc_total_samples_used > perc_data_to_use)\
             or (perc_incl_data_reused is not None and\
             perc_included_samples_reused > perc_incl_data_reused):
            upper_bound = est_n_cells_per_metacell
        elif (perc_data_to_use is not None and\
             perc_total_samples_used < perc_data_to_use)\
             or (perc_incl_data_reused is not None and\
             perc_included_samples_reused < perc_incl_data_reused):
            lower_bound = est_n_cells_per_metacell
        else:
            raise ValueError("Error in optimize_knn_by_size_with_perc_data")

    n_cells_per_metacell = est_n_cells_per_metacell
    sample_indices = temp_sample_indices
    return n_cells_per_metacell, sample_indices

def get_sample_indices_by_optimizing_knn_by_size_with_perc_data(
    adata,
    pca_slot,
    dist_slot,
    size,
    perc_data_to_use,
    perc_incl_data_reused,
    exact_size,
    seed,
    verbose,
    njobs
):
    pca_array = adata.obsm[pca_slot].copy()

    if dist_slot not in list(adata.obsp.keys()):
        warnings.warn(dist_slot + " not in adata.obsp. Computing correlation distance...")
        corr_distance(adata,
                      use_reduction=True,
                      reduction_slot=pca_slot,
                      key_added='corr_dist')
        dist_slot = 'corr_dist'

    knn_groups_indices_df = get_mc_indices_by_condensing(pca_array, size, exact_size, seed, verbose, njobs)
    neighbors_knn(adata,
                          max_knn=int(np.ceil(adata.shape[0]/np.sqrt(size))),#Conservative upper bound on KNN, #None
                          dist_slot=dist_slot,
                          key_added='knn4MC',
                          njobs = 1)
    knn_array = adata.uns['knn4MC']

    n_cells_per_metacell, temp_sample_indices = optimize_knn_by_size_with_perc_data(
        knn_groups_indices_df,
        knn_array,
        perc_data_to_use,
        perc_incl_data_reused,
        verbose
    )
    adata.obs["selected_for_mc"] = pd.Categorical(np.isin(np.arange(adata.shape[0]), temp_sample_indices))
    return adata, size, n_cells_per_metacell

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ---------------------------------------------------------------------------
# ----------------- get_sample_indices_by_optimizing_params -----------------
# ---------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def get_sample_indices_by_optimizing_params(
    adata,
    pca_slot,
    dist_slot,
    size,
    n_cells_per_metacell,
    perc_data_to_use,
    perc_incl_data_reused,
    exact_size,
    seed,
    verbose,
    njobs
):
    if perc_data_to_use is not None and perc_incl_data_reused is not None:
        raise ValueError("perc_data_to_use and perc_incl_data_reused contradict: one must be set to None")

    if size is None:
        return get_sample_indices_by_optimizing_size_by_knn_with_perc_data(
            adata,
            pca_slot,
            dist_slot,
            n_cells_per_metacell,
            perc_data_to_use,
            perc_incl_data_reused,
            exact_size,
            seed,
            verbose,
            njobs
        )
    elif n_cells_per_metacell is None:
        return get_sample_indices_by_optimizing_knn_by_size_with_perc_data(
            adata,
            pca_slot,
            dist_slot,
            size,
            perc_data_to_use,
            perc_incl_data_reused,
            exact_size,
            seed,
            verbose,
            njobs
        )
    else:
        raise ValueError("size and n_cells_per_metacell contradict: one must be set to None")

    adata.obs["selected_for_mc"] = pd.Categorical(np.isin(np.arange(adata.shape[0]), sample_indices))
    return adata, size, n_cells_per_metacell

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ---------------------------------------------------------------------------
# ------------------------------ MAIN FUNCTION ------------------------------
# ---------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def _representative_metacells(
    adata,
    counts = None,
    pca_slot = "X_pca",
    dist_slot = "corr_dist",
    size = 250,
    n_cells_per_metacell = None,
    min_median_depth = 10000,
    perc_data_to_use = None,
    perc_incl_data_reused = None,
    exact_size = True,
    seed = 0,
    key_added = "metacells",
    verbose = True,
    njobs = 1,
    smart_sample = True
):
    # Like a ROC plot as you increase % total samples used you also start to increase
    # % sample reuse. It would be nice if we could somehow get the curve and then users
    # can specify what point on the curve they want. Do they want to be strict about reusing
    # samples or do they prefer using as many samples as possible despite sample reuse? Or
    # do they just want to set a minimum median depth and the point will fall wherever?
    # It would be great if they could have that option.
    counts = get_counts_as_df(counts, adata)

    if adata.shape[0] < size:
        raise ValueError("Number of metacells requested (" + str(size) + \
                         ") is greater than the number of samples in adata (" +  \
                         str(adata.shape[0]) + ").")

    if np.sum([size != None,
               min_median_depth != None,
               n_cells_per_metacell != None,
               perc_data_to_use != None,
               perc_incl_data_reused != None]) != 2:
        print("Currently set parameters:")
        if size is not None: print("\tsize=" + str(size))
        if min_median_depth is not None: print("\tmin_median_depth=" + str(min_median_depth))
        if n_cells_per_metacell is not None: print("\tn_cells_per_metacell=" + str(n_cells_per_metacell))
        if perc_data_to_use is not None: print("\tperc_data_to_use=" + str(perc_data_to_use))
        if perc_incl_data_reused is not None: print("\tperc_incl_data_reused=" + str(perc_incl_data_reused))

        raise ValueError("Exactly two of the following options must be set:" + \
                         "\n\tsize, min_median_depth, n_cells_per_metacell, perc_data_to_use, perc_incl_data_reused." + \
                         "\n\tNote: min_median_depth & n_cells_per_metacell cannot both be set." + \
                         "\n\tNote: perc_data_to_use & perc_incl_data_reused cannot both be set.")
    if np.sum([min_median_depth != None, n_cells_per_metacell != None]) == 2:
        raise ValueError("One or none but not both of the following options can be set:" + \
                         "\nmin_median_depth and n_cells_per_metacell. " + \
                         "Leave the other option as None.")
    if np.sum([perc_data_to_use != None, perc_incl_data_reused != None]) == 2:
        raise ValueError("One or none but not both of the following options can be set:" + \
                         "\nperc_data_to_use and perc_incl_data_reused. " + \
                         "Leave the other option as None.")
    # ADD CODE TO DISPLAY WHAT IS CURRENTLY SET:

    if n_cells_per_metacell is None:
        if min_median_depth is not None:
            ave_cell_depth = np.mean(np.sum(counts.values, axis = 1))
            n_cells_per_metacell = int(np.ceil(min_median_depth / ave_cell_depth))
    else:
        n_cells_per_metacell = int(n_cells_per_metacell)

    if perc_data_to_use is None and perc_incl_data_reused is None:
        adata, size, n_cells_per_metacell = get_sample_indices_with_preset_params(
            adata,
            counts,
            pca_slot,
            dist_slot,
            size,
            n_cells_per_metacell,
            min_median_depth,
            exact_size,
            seed,
            verbose,
            njobs,
            smart_sample
        )
    else:
        adata, size, n_cells_per_metacell = get_sample_indices_by_optimizing_params(
            adata,
            pca_slot,
            dist_slot,
            size,
            n_cells_per_metacell,
            perc_data_to_use,
            perc_incl_data_reused,
            exact_size,
            seed,
            verbose,
            njobs
        )

    if verbose:
        print("Finalized parameters:")
        print("\tsize=" + str(size))
        print("\tmin_median_depth=" + str(min_median_depth))
        print("\tn_cells_per_metacell=" + str(n_cells_per_metacell))

    knn_array = adata.uns['knn4MC']
    sample_indices = np.where(adata.obs["selected_for_mc"].values == True)[0]

    if verbose and pca_slot in list(adata.obsm.keys()):
        sc.pl.pca(adata, color = "selected_for_mc")

    # Use the KNN to identify the neighbors of each selected index
    metacell_knn = knn_array[sample_indices,:n_cells_per_metacell].copy()
    metacell_array = np.sum(counts.values[metacell_knn,:], axis=1)
    metacell_df = pd.DataFrame(metacell_array,
                               columns = counts.columns,
                               index = ["metacell_" + str(i) for i in range(len(sample_indices))])

    # Output statistics
    # Percent total of data used
    prop_total_samples_used = len(np.unique(metacell_knn.flatten()))/adata.shape[0]
    perc_total_samples_used = np.round(prop_total_samples_used*100,1)
    if verbose: print(str(perc_total_samples_used) + "% total samples were used.")
    # Percent total of data used more than once
    n_times_sample_used = np.unique(metacell_knn.flatten(), return_counts = True)[1]
    prop_samples_used_mult_times = np.mean(n_times_sample_used > 1)
    perc_samples_used_mult_times = np.round(prop_samples_used_mult_times*100,1)
    if verbose: print(str(perc_samples_used_mult_times) + "% of these included samples were used more than once.")
    mean_n_times_used = np.round(np.mean(n_times_sample_used),1)
    stdev_n_times_used = np.round(np.std(n_times_sample_used),1)
    if verbose: print("Sample were used an average of " + str(mean_n_times_used) + " times (stdev=" + str(stdev_n_times_used) + ")")

    mean_n_times_used = np.round(np.mean(n_times_sample_used[n_times_sample_used > 1]),1)
    med_n_times_used = np.round(np.median(n_times_sample_used[n_times_sample_used > 1]),1)
    stdev_n_times_used = np.round(np.std(n_times_sample_used[n_times_sample_used > 1]),1)
    if verbose: print("Reused samples samples were used an average of " + str(mean_n_times_used) + " times (stdev=" + str(stdev_n_times_used) + ")")

    depth_values = np.sum(metacell_df, axis = 1).values
    if verbose:
        print("Depth Statistics")
        print("\tMean: " + str(np.round(np.mean(depth_values),1)))
        print("\tMedian: " + str(np.round(np.median(depth_values),1)))
        print("\tStdev: " + str(np.round(np.std(depth_values),1)))
        print("\tMin: " + str(np.min(depth_values)))
        print("\tMax: " + str(np.max(depth_values)))

    metacell_df.attrs['size'] = size
    metacell_df.attrs['n_cells_per_metacell'] = n_cells_per_metacell
    metacell_df.attrs['min_median_depth'] = min_median_depth
    metacell_df.attrs['median_depth'] = np.median(depth_values)
    metacell_df.attrs['perc_total_samples_used'] = prop_total_samples_used*100
    metacell_df.attrs['perc_included_samples_used_mult_times'] = prop_samples_used_mult_times*100
    metacell_df.attrs['mean_n_times_samples_used'] = np.mean(n_times_sample_used)
    metacell_df.attrs['stdev_n_times_samples_used'] = np.std(n_times_sample_used)

    adata.uns[key_added] = metacell_df

    # Output plots
    if verbose:
        metacell_anndata = anndata.AnnData(metacell_df)
        metacell_anndata.var["mt"] = metacell_anndata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            metacell_anndata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        sc.pl.violin(
            metacell_anndata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
        )

def _representative_metacells_multiclusters(
    adata,
    counts = None,
    pca_slot = "X_pca",
    dist_slot = "corr_dist",
    clusters_slot=None,
    size = 250,
    n_cells_per_metacell = None,
    min_median_depth = 10000,
    perc_data_to_use = None,
    perc_incl_data_reused = None,
    exact_size = True,
    seed = 0,
    key_added = "metacells",
    verbose = True,
    njobs = 1,
    smart_sample = True,
    copy = False
):
    if copy: adata = adata.copy()
    counts = get_counts_as_df(counts, adata)

    if clusters_slot is not None:
        unique_clusters = np.unique(adata.obs[clusters_slot].values)
        n_unique_clusters = len(unique_clusters)
        for i in tqdm(range(n_unique_clusters), desc="cluster metacells") if verbose else range(n_unique_clusters):
            clust = unique_clusters[i]
            clust_size = np.sum(adata.obs[clusters_slot]==clust)
            if clust_size < size:
                warnings.warn("Number of samples (" + str(clust_size) + ") in cluster " + \
                              str(clust) + " is less than number of metacells requested (" + \
                              str(size) + "). Skipping cluster.")
                continue

            adata_clust = adata[adata.obs[clusters_slot]==clust,].copy()
            counts_clust = counts[adata.obs[clusters_slot]==clust].copy()

            _representative_metacells(
                adata_clust,
                counts_clust,
                pca_slot,
                dist_slot,
                size,
                n_cells_per_metacell,
                min_median_depth,
                perc_data_to_use,
                perc_incl_data_reused,
                exact_size,
                seed,
                key_added,
                verbose = False,
                njobs = njobs,
                smart_sample = smart_sample
            )
            clust_metacells = adata_clust.uns[key_added]

            adata.uns[str(key_added + "_" + str(clust))] = clust_metacells
    else:
        _representative_metacells(
            adata,
            counts,
            pca_slot,
            dist_slot,
            size,
            n_cells_per_metacell,
            min_median_depth,
            perc_data_to_use,
            perc_incl_data_reused,
            exact_size,
            seed,
            key_added,
            verbose,
            njobs,
            smart_sample
        )

    if copy: return adata
