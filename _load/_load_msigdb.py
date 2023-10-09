### Import dependencies
from ._load_get_path import __get_pyviper_dir
from ..interactome import Interactome

### EXPORT LIST
__all__ = ['load_msigdb_regulon']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ----------------------- LOAD MSIGDB REGULONS FUNCTIONS ----------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------

def __merge_dicts(dict1, dict2):
    res = {**dict1, **dict2}
    return res
def __get_msigdb_reg_path(collection):
    reg_path = __get_pyviper_dir() + "/data/regulons/msigdb-" + collection + "-as-regulon.tsv"
    return(reg_path)
def __load_msigdb_from_tsv(collection):
    reg_path = __get_msigdb_reg_path(collection.lower())
    reg = Interactome("MSigDB_" + collection.upper(), reg_path)
    return(reg)
def __load_msigdb_regulon_single(collection = "c2"):
    reg = None
    if(collection.lower() in ["c2", "c5", "c6", "c7", "h"]):
        reg = __load_msigdb_from_tsv(collection)
    return(reg)
def __load_msigdb_regulon_multiple(collection = ["h", "c2"]):
    combined_dict = {}
    for i in range(len(collection)):
        new_dict = __load_msigdb_regulon_single(collection[i]).regDict
        combined_dict = __merge_dicts(combined_dict, new_dict)
    combined_interactome = Interactome(name = '+'.join(collection),
                                       regDict = combined_dict)
    return(combined_interactome)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------

def load_msigdb_regulon(collection = "h"):
    reg = None
    if(type(collection) is str):
        reg = __load_msigdb_regulon_single(collection)
    elif(type(collection) is list):
        reg = __load_msigdb_regulon_multiple(collection)
    return(reg)
