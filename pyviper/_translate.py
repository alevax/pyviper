### ---------- IMPORT DEPENDENCIES ----------
from tqdm import tqdm
from ._load._load_translate import load_human2mouse
import numpy as np
import pandas as pd
from anndata import AnnData
from warnings import warn

### ---------- EXPORT LIST ----------
__all__ = ['translate_adata_index', '_detect_name_type']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def _detect_name_type(input_array):
    nrow = len(input_array)
    if nrow == 0: return(None) #redundant: case handled by 0-length for loop
    found_match = False
    human2mouse = load_human2mouse()
    gene_name_format = None
    i = 0
    for i in range(nrow):
        gene = str(input_array[i])
        if(gene in human2mouse["mouse_symbol"].values):
            gene_name_format = "mouse_symbol"
            break
        elif(gene in human2mouse["human_symbol"].values):
            gene_name_format = "human_symbol"
            break
        elif(gene in human2mouse["mouse_ensembl"].values):
            gene_name_format = "mouse_ensembl"
            break
        elif(gene in human2mouse["human_ensembl"].values):
            gene_name_format = "human_ensembl"
            break
        elif(gene in human2mouse["mouse_entrez"].values):
            gene_name_format = "mouse_entrez"
            break
        elif(gene in human2mouse["human_entrez"].values):
            gene_name_format = "human_entrez"
            break
    return(gene_name_format)

def _translate_genes_array(current_gene_names, desired_format):
    # if desired_format in ['human_symbol', 'human_ensembl', 'human_entrez']:
    #     translate_df = load_mouse2human()
    # elif desired_format in ['mouse_symbol', 'mouse_ensembl', 'mouse_entrez']:
    #     translate_df = load_human2mouse()
    if desired_format in ['human_symbol', 'human_ensembl', 'human_entrez',
                          'mouse_symbol', 'mouse_ensembl', 'mouse_entrez']:
        translate_df = load_human2mouse()
    else:
        raise ValueError("Error: desired_format is not one the following:"
                         + "\n\t\t mouse_symbol, mouse_ensembl, mouse_entrez,"
                         + "\n\t\t human_symbol, human_ensembl, human_entrez")

    current_format = _detect_name_type(current_gene_names)
    if current_format is None:
        raise ValueError("Error: could not detect current_format as one of the following:"
                         + "\n\t\t mouse_symbol, mouse_ensembl, mouse_entrez,"
                         + "\n\t\t human_symbol, human_ensembl, human_entrez")

    translate_df = translate_df.sort_values(by=current_format, ascending=True)
    dict_column_current_format = translate_df[current_format].values
    dict_column_desired_format = translate_df[desired_format].values

    # Get positions of translated values
    positions = np.searchsorted(dict_column_current_format, current_gene_names)

    # Identify input values that have no translation available
    # These two lines of code are essentially:
        # mask = (positions < len(dict_column_current_format)) & (dict_column_current_format[positions] == current_gene_names)
    # But dict_column_current_format[positions] causes error when positions = len(dict_column_current_format), which
    # happens if we have elements sorted to the end, such as with "A", "B", "C" with element "D".
    mask = (positions < len(dict_column_current_format))
    selected_elements = mask[mask]
    selected_elements[dict_column_current_format[positions[mask]] != current_gene_names[mask]] = False
    if desired_format == 'mouse_entrez':
        selected_elements[dict_column_desired_format[positions[mask]] == -1] = False
    if current_format == 'mouse_entrez':
        selected_elements[dict_column_current_format[positions[mask]] == -1] = False
    mask[mask] = selected_elements

    # Get translated values
    translation = np.array([None] * len(current_gene_names))
    translation[mask] = dict_column_desired_format[positions[mask]]

    return translation

def keep_first_duplicate_strings(arr):
    arr[arr == None] = '-1'
    _, index = np.unique(arr, return_index=True)
    result = np.full_like(arr, fill_value=np.nan, dtype=object)
    result[index] = arr[index]
    result[result == '-1'] = np.nan
    return result

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def translate_adata_index(adata, desired_format, eliminate = True, copy = False):
    """\
    Take adata.var.index.names, replace them with a translation of desired_format,
    and move the original names to a new column in var. The current name format
    and desired_format of the gene names in adata.var.index should be one of the
    following:
        mouse_symbol, mouse_ensembl, mouse_entrez, human_symbol, human_ensembl, or human_entrez

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata object.
    desired_format
        Desired format can be one of six strings: "mouse_symbol", "mouse_ensembl",
        "mouse_entrez", "human_symbol", "human_ensembl", or "human_entrez".
    eliminate (default: True)
        Whether to eliminate var rows that don't have a translation. Otherwise,
        None will be left in place in the index.
    copy (default: False)
        Whether to return a translated copy (True) or to instead translate the
        original input (False).

    Returns
    -------
    The translated anndata object.
    """
    # So all translation happen on a new adata that is returned. Otherwise, we
    # edit the original, but eliminate doesn't work if the output isn't taken.
    if copy: adata = adata.copy()
    current_format = _detect_name_type(adata.var.index.values)
    adata.var[current_format] = adata.var.index.values.astype(str)
    adata.var[desired_format] = _translate_genes_array(adata.var[current_format], desired_format)
    # return adata
    adata.var[desired_format] = keep_first_duplicate_strings(adata.var[desired_format].values)
    if eliminate:
        adata._inplace_subset_var(~pd.isna(adata.var[desired_format]))
        adata._inplace_subset_var(adata.var[desired_format] != "nan")
        adata._inplace_subset_var(adata.var[desired_format] != "NaN")

    adata.var.set_index(desired_format, inplace=True)

    if copy: return adata


def translate(adata, desired_format, eliminate = True, copy = False):
    """\
    Performs translation of an AnnData, np.ndarray, or list object. For AnnData,
    takes adata.var.index.names, replaces them with a translation of desired_format,
    and moves the original names to a new column in var. The current name format
    and desired_format of the gene names in adata.var.index should be one of the
    following:
        mouse_symbol, mouse_ensembl, mouse_entrez, human_symbol, human_ensembl, or human_entrez

    Parameters
    ----------
    adata
        Gene expression, protein activity or pathways stored in an anndata
        object. Alternate inputs: np.ndarray or list.
    desired_format
        Desired format can be one of six strings: "mouse_symbol", "mouse_ensembl",
        "mouse_entrez", "human_symbol", "human_ensembl", or "human_entrez".
    eliminate (default: True)
        Whether to eliminate var rows that don't have a translation. Otherwise,
        None will be left in place in the index.
    copy (default: False)
        Whether to return a translated copy (True) or to instead translate the
        original AnnData (False). True when given np.ndarray or list.

    Returns
    -------
    Returns a translated object when copy = True or input is np.ndarray or list.
    Modifies the original when copy = False and input is AnnData.
    """
    if isinstance(adata, AnnData):
        if copy is True:
            return translate_adata_index(adata, desired_format, eliminate, copy)
        else:
            translate_adata_index(adata, desired_format, eliminate, copy)
    elif isinstance(adata, np.ndarray):
        if copy is False:
            warn("np.ndarray supplied:" +
                 " overriding copy=False and returning a translated copy.")
        original_shape = adata.shape
        adata = adata.flatten()
        current_format = _detect_name_type(adata)
        adata = _translate_genes_array(adata, desired_format)
        # Switch to "-1" so shape isn't affected
        adata[adata == None] = '-1'
        adata[adata == np.nan] = '-1'
        adata = keep_first_duplicate_strings(adata)
        adata = adata.reshape(original_shape)
        if eliminate:
            adata = adata[~pd.isna(adata)]
            adata = adata[adata != "nan"]
            adata = adata[adata != "NaN"]
        return adata
    elif isinstance(adata, list):
        if copy is False:
            warn("list supplied to translate:" +
                 " overriding copy=False and returning a translated copy.")
        adata = np.array(adata)
        current_format = _detect_name_type(adata)
        adata = _translate_genes_array(adata, desired_format)
        adata = keep_first_duplicate_strings(adata)
        if eliminate:
            adata = adata[~pd.isna(adata)]
            adata = adata[adata != "nan"]
            adata = adata[adata != "NaN"]
        adata = list(adata)
        return adata
