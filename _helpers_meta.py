### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np


def get_all_regs(mat_list):
    # Iterate through each DataFrame to collect unique gene names
    all_regs_set = set()
    for mat in mat_list:
        reg_names = mat.columns.tolist()
        all_regs_set.update(reg_names)

    # Sort the gene names alphabetically
    all_regs = sorted(list(all_regs_set))

    return all_regs

def get_resized_mats(mat_list, empty_value = np.nan):
    # Get names of all regulators
    all_regs = get_all_regs(mat_list)

    # Iterate through each DataFrame to normalize it
    resized_mat_list = []
    for mat in mat_list:
        # Create a DataFrame with empty values (e.g. NaN or 0) for missing gene names
        mat = mat.copy()
        missing_regs = list(set(all_regs) - set(mat.columns))
        empty_df = pd.DataFrame(empty_value, index=mat.index, columns=list(missing_regs))
        # Concatenate the original DataFrame with the empty (e.g. NaN or 0) DataFrame
        resized_mat = pd.concat([mat, empty_df], axis=1)
        # Sort the columns alphabetically
        resized_mat = resized_mat[all_regs]
        resized_mat_list.append(resized_mat)

    return resized_mat_list
