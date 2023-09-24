### Import dependencies
import pandas as pd
import pathlib
from .interactome import Interactome

### EXPORT LIST
__all__ = []

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ GENERAL FUNCTION -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def __get_pyther_dir():
    pyther_dir = str(pathlib.Path(__file__).parent.resolve())
    return(pyther_dir)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------ LOAD TRANSLATE DFS FUNCTIONS -----------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def load_mouse2human():
    mouse2human = pd.read_csv(__get_pyther_dir() + "/data/translate/human2mouse.csv")
    del mouse2human[mouse2human.columns[0]]
    return(mouse2human)
def load_human2mouse():
    human2mouse = pd.read_csv(__get_pyther_dir() + "/data/translate/mouse2human.csv")
    del human2mouse[human2mouse.columns[0]]
    return(human2mouse)

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# --------------------------- LOAD REGULATORS FUNCTIONS --------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------

def __load_regulators(path_to_txt):
    with open(path_to_txt) as temp_file:
        regulator_set = [line.rstrip('\n') for line in temp_file]
    return(regulator_set)

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

def load_TFs(path_to_tfs = None, species = "human"):
    if path_to_tfs is None:
        if species == "human":
            path_to_tfs = __get_pyther_dir() + "/data/regulatorIDs/tfs-hugo.txt"
        elif species == "mouse":
            path_to_tfs = __get_pyther_dir() + "/data/regulatorIDs/tfs-hugo-mouse.txt"
    tfs_list = __load_regulators(path_to_tfs)
    return(tfs_list)

def load_coTFs(path_to_cotfs = None, species = "human"):
    if path_to_cotfs is None:
        if species == "human":
            path_to_cotfs = __get_pyther_dir() + "/data/regulatorIDs/cotfs-hugo.txt"
        elif species == "mouse":
            path_to_cotfs = __get_pyther_dir() + "/data/regulatorIDs/cotfs-hugo-mouse.txt"
    cotfs_list = __load_regulators(path_to_cotfs)
    return(cotfs_list)

def load_sig(path_to_sig = None, species = "human"):
    if path_to_sig is None:
        if species == "human":
            path_to_sig = __get_pyther_dir() + "/data/regulatorIDs/sig-hugo.txt"
        elif species == "mouse":
            path_to_sig = __get_pyther_dir() + "/data/regulatorIDs/sig-hugo-mouse.txt"
    sig_list = __load_regulators(path_to_sig)
    return(sig_list)

def load_surf(path_to_surf = None, species = "human"):
    if path_to_surf is None:
        if species == "human":
            path_to_surf = __get_pyther_dir() + "/data/regulatorIDs/surface-hugo.txt"
        elif species == "mouse":
            path_to_surf = __get_pyther_dir() + "/data/regulatorIDs/surface-hugo-mouse.txt"
    surf_list = __load_regulators(path_to_surf)
    return(surf_list)

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
    reg_path = __get_pyther_dir() + "/data/regulons/msigdb-" + collection + "-as-regulon.tsv"
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
