### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import anndata

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def mat_to_anndata(mat):
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




# # Function assumes features as rows and observations as columns
# # Numerator Functions:
#     # Median - numpy.median
#     # Mean - numpy.mean
# # Denominator Functions:
#     # Median absolute deviation - mad_from_R
#     # Standard deviation - statistics.stdev
# def rank_norm(x, NUM_FUN=np.median, DEM_FUN = mad_from_R, trim=0):
#     rank = rankdata(x, axis=0)
#     median = NUM_FUN(rank, axis=1, keepdims=True)#np.median(rank, axis=1, keepdims=True)
#     mad = np.apply_along_axis(DEM_FUN, 1, rank)
#
#     x = ((rank - median)/mad[:, np.newaxis])
#
#     print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
#     print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
#     print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))
#     print("- Features to Remove:")
#
#     # Take care of infinite values
#     max_finite = np.nanmax(x[np.isfinite(x)])
#     min_finite = np.nanmin(x[np.isfinite(x)])
#     x[np.isposinf(x)] = max_finite
#     x[np.isneginf(x)] = min_finite
#
#     x = np.where(np.isnan(x), np.nanmin(x), x)
#     x = np.clip(x, a_min=np.nanmin(x), a_max=np.nanmax(x))
#     print("- Removing NULL/NA features ...")
#     x = x[~np.isnan(x).any(axis=1)]
#
#     print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
#     print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
#     print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))
#
#     return x
# def mad_from_R(x, center=None, constant=1.4826, low=False, high=False):
#     if center is None:
#         center=np.median(x)
#     x = x[~np.isnan(x)] if np.isnan(x).any() else x
#     n = len(x)
#     if (low or high) and n % 2 == 0:
#         if low and high:
#             raise ValueError("'low' and 'high' cannot be both True")
#         n2 = n // 2 + int(high)
#         return constant * np.sort(np.abs(x - center))[n2]
#     return constant * np.median(np.abs(x - center))




# def slice_concat(inner_function, gex_data ,bins = 10, write_local = True, **kwargs):
#     #kwargs are the parameters for the inner function.
#     #slice the data cells * genes
#
#     result_list = []
#     size = int(gex_data.shape[0]/bins)
#     residue = gex_data.shape[0] % bins
#
#     if write_local:
#         os.mkdir('temp')
#
#         for i in range(bins-1):
#             segment = gex_data[i*size: i*size + size,]
#             temp_result = inner_function(segment, **kwargs)
#
#             if type(temp_result) == anndata._core.anndata.AnnData:
#                 temp_result = temp_result.to_df()
#
#
#             temp_result.to_csv('temp/'+ str(i) + '.csv')
#
#         # the last one
#         segment = gex_data[(bins-1)*size: bins*size + residue,]
#         temp_result = inner_function(segment, **kwargs)
#
#         if type(temp_result) == anndata._core.anndata.AnnData:
#             temp_result = temp_result.to_df()
#
#         temp_result.to_csv('temp/'+ str(bins-1) + '.csv')
#
#
#         all_file_list=os.listdir('temp')
#         for single_file in all_file_list:
#             result_list.append(pd.read_csv(os.path.join('temp',single_file)))
#
#         shutil.rmtree('temp')
#
#     else:
#         for i in range(bins):
#             segment = gex_data[i*size: i*size + size,]
#             result_list.append(inner_function(segment, **kwargs))
#
#
#     # concat result
#
#     result = pd.concat(result_list,axis=0).reset_index(drop = True)
#     result.set_index(keys='index',inplace=True)
#     return result
