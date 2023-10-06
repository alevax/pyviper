### ---------- IMPORT DEPENDENCIES ----------
from tqdm import tqdm

### ---------- EXPORT LIST ----------
__all__ = ['translate_adata_index']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# ---------------------------- TRANSLATE SMALL PICTURE FUNCS ---------------------------
def translate_gene_name(input_gene_name,
              translate_df,
              current_format = "mouse_symbol",
              desired_format = "human_symbol"):
    translate_df_input_gene_row = translate_df.loc[translate_df[current_format] == input_gene_name]
    if(translate_df_input_gene_row.empty):
        output_gene_name = None
    else:
        output_gene_name = translate_df_input_gene_row[desired_format].values[0]
    return(output_gene_name)
def translate_adata_index_into_new_var_column(
         adata,
         translate_df,
         current_format = "mouse_symbol",
         desired_format = "human_symbol",
         show_progress_bar = True):
    current_gene_names = adata.var.index.values
    n_genes = len(current_gene_names)
    desired_gene_names = [None]*n_genes

    iterator = (i for i in range(n_genes)) if not show_progress_bar else tqdm(range(n_genes))
    for i in iterator:
        desired_gene_names[i] = translate_gene_name(
            current_gene_names[i],
            translate_df,
            current_format,
            desired_format
        )
    adata.var[desired_format] = desired_gene_names
    return(adata)
# ---------------------------- TRANSLATE BIG PICTURE FUNCS ---------------------------
def translate_adata_index_from_to_with_translate_df(adata,
         translate_df,
         current_format = "mouse_symbol",
         desired_format = "human_symbol",
         show_progress_bar = True):
    adata = translate_adata_index_into_new_var_column(
         adata,
         translate_df,
         current_format,
         desired_format,
         show_progress_bar)
    adata.var[current_format] = adata.var.index.values
    adata.var.set_index(desired_format, inplace=True)
    return(adata)
def translate_adata_index_from_to(adata,
                                  current_format = "mouse_symbol",
                                  desired_format = "human_symbol",
                                  show_progress_bar = True):
    acceptable_formats = ["mouse_symbol", "mouse_ensembl", "human_symbol", "human_ensembl"]
    if(current_format in ["mouse_symbol", "mouse_ensembl"]):
        translate_df = load_mouse2human()
    elif(current_format in ["human_symbol", "human_ensembl"]):
        translate_df = load_human2mouse()
    if current_format not in acceptable_formats:
        raise ValueError("Error: index of adata.var is not one the following:"
                         + "\n\t\t mouse_symbol, mouse_ensembl, human_symbol, human_ensembl")
    if(desired_format in adata.var.columns):
        adata.var[current_format] = adata.var.index
        adata.var.set_index(desired_format, inplace=True)
    else:
        adata = translate_adata_index_from_to_with_translate_df(adata,
             translate_df,
             current_format,
             desired_format = desired_format,
             show_progress_bar = show_progress_bar)
    return(adata)


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def translate_adata_index(adata,
                          desired_format = "human_symbol",
                          current_format = None,
                          show_progress_bar = True):
    acceptable_formats = ["mouse_symbol", "mouse_ensembl", "human_symbol", "human_ensembl"]
    if desired_format not in acceptable_formats:
        raise ValueError("Error: desired_format is not one the following:"
                         + "\n\t\t mouse_symbol, mouse_ensembl, human_symbol, human_ensembl")
    if current_format is None:
        current_format = detect_index_name_type(adata)
    adata = translate_adata_index_from_to(adata,
                                          current_format,
                                          desired_format,
                                          show_progress_bar)
    return(adata)
