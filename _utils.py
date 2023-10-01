### ---------- IMPORT DEPENDENCIES ----------
import anndata
import pandas as pd
import numpy as np
from ._helpers import *
from .aREA.aREA_classic import *
from ._translate import *
from .load import load_human2mouse, load_mouse2human
# from .interactome import Interactome

### ---------- EXPORT LIST ----------
__all__ = ["path_enr", "compute_cluster_stouffer_anndata", "compute_cluster_stouffer"]

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# ---------------------------- DETECT FORMAT TYPE FUNCS ---------------------------
def __detect_interactome_name_type(interactome):
    return(__detect_name_type(np.array(list(interactome.get_targetSet()))))
def __detect_index_name_type(adata):
    return(__detect_name_type(adata.var.index))
def __detect_name_type(input_array):
    nrow = len(input_array)
    if nrow == 0: return(None) #redundant: case handled by 0-length for loop
    found_match = False
    human2mouse = load_human2mouse()
    gene_name_format = None
    i = 0
    for i in range(nrow):
        gene = input_array[i]
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
    return(gene_name_format)


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def path_enr(adata,
             interactome,
             layer = None,
             verbose = True,
             transfer_obs = True):
    """\
    Allows the individual to infer normalized enrichment scores of pathways
    using the analytical ranked enrichment analysis (aREA) function.

    This is a wrapper for the aREA function in Pyther.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    interactome
        The interactome object containing pathways.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
    """
    if(verbose): print("Checking interactome names format...")
    interactome_gene_name_format = __detect_interactome_name_type(interactome)
    if(verbose): print("Interactome names are formatted as " + interactome_gene_name_format + ".")
    if(verbose): print("Checking adata names format...")
    adata_gene_name_format_original = __detect_index_name_type(adata)
    if(verbose): print("adata names are formatted as " + adata_gene_name_format_original + ".")
    make_adata_names_format_match_interactome = False

    if(interactome_gene_name_format != adata_gene_name_format_original):
        make_adata_names_format_match_interactome = True
        if(verbose): print("Translating adata names to match interactome...")
        adata = translate_adata_index(adata,
                                      current_format = adata_gene_name_format_original,
                                      desired_format = interactome_gene_name_format)

    # aREA takes the pathways interactome and the adata
    if(verbose): print("Running aREA using to calculate pathway enrichment...")
    interactome.filter_targets(adata.var_names)
    path_enr_mat = aREA_classic(adata, interactome, eset_filter = False, layer = layer, min_targets=0, verbose = verbose)

    if(make_adata_names_format_match_interactome is True):
        if(verbose): print("Returning adata names to original state...")
        adata = translate_adata_index(adata,
                                      current_format = interactome_gene_name_format,
                                      desired_format = adata_gene_name_format_original)
    # Create a new Anndata object
    pwe_data = mat_to_anndata(path_enr_mat)
    # This means we did pathway enrichment on VIPER: adata is pax_data
    if hasattr(adata, "gex_data"):
        pwe_data.gex_data = adata.gex_data
        adata.gex_data = None
        pwe_data.pax_data = adata
    # This means we did pathway enrichment on gex: adata is gex_data
    else:
        pwe_data.gex_data = adata

    if transfer_obs is True:
        pwe_data.obs = adata.obs

    return(pwe_data)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ----------------------------- STOUFFER FUNCTIONS ----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def compute_cluster_stouffer_anndata(adata, obs_column_name, layer = None):
    if layer is None:
        dat_df = pd.DataFrame(adata.X,
                              index=adata.obs_names,
                              columns=adata.var_names)
    else:
        dat_df = pd.DataFrame(adata.layers[layer],
                              index=adata.obs_names,
                              columns=adata.var_names)
    cluster_vector = adata.obs[obs_column_name]
    result_df = compute_cluster_stouffer(dat_df, cluster_vector)
    return anndata.AnnData(result_df)
    # return mat_to_anndata(result_df)

def compute_cluster_stouffer(dat_df, cluster_vector):
    # Ensure cluster_vector has the same number of samples as rows in dat_df
    if len(cluster_vector) != dat_df.shape[0]:
        raise ValueError("Cluster vector length does not match the number of rows in the DataFrame.")

    # Convert the DataFrame to a NumPy array
    dat_array = dat_df.to_numpy()

    # Find unique clusters and initialize arrays to store Stouffer scores
    unique_clusters, cluster_indices = np.unique(cluster_vector, return_inverse=True)
    n_clusters = len(unique_clusters)
    n_genes = dat_df.shape[1]
    stouffer_scores = np.zeros((n_clusters, n_genes))

    # Calculate the denominator for Stouffer scores for each cluster
    cluster_sizes = np.bincount(cluster_indices)
    sqrt_cluster_sizes = np.sqrt(cluster_sizes)

    # Calculate Stouffer scores for each cluster and gene
    for i in range(n_clusters):
        cluster_mask = (cluster_indices == i)
        cluster_data = dat_array[cluster_mask]
        stouffer_scores[i, :] = np.sum(cluster_data, axis=0) / sqrt_cluster_sizes[i]

    # Create a DataFrame from the computed Stouffer scores
    result_df = pd.DataFrame(stouffer_scores, index=unique_clusters, columns=dat_df.columns)

    return result_df
