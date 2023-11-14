### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import anndata
from statsmodels.stats import multitest

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def mat_to_anndata(mat):
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
                               var=mat_features,
                               dtype=np.float64)
    return(pax_data)

def _adjust_p_values(p_values):
	# Helper function for *nes_to_pval* in module "tl" 
    # correct p values with FDR
    _, corrected_p_values, _, _ = multitest.multipletests(p_values, method='fdr_bh')
    return corrected_p_values
    




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
