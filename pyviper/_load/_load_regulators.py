### Import dependencies
from ._load_get_path import __get_pyviper_dir
from ..config import config

### EXPORT LIST
__all__ = ['load_TFs', 'load_coTFs', 'load_sig', 'load_surf']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# --------------------------- LOAD REGULATORS FUNCTIONS --------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------

def __load_regulators_from_txt(path_to_txt):
    with open(path_to_txt) as temp_file:
        regulator_set = [line.rstrip('\n') for line in temp_file]
    return(regulator_set)

def _load_regulators(group, species = None, path_to_regs = None):
    if species is None:
        species = config['regulators_species']
    if path_to_regs is None:
        config_path_to_regs = config['regulators_filepaths'][species][group]
        if config_path_to_regs is None:
            if species == "human":
                config_path_to_regs = __get_pyviper_dir() + "/data/regulatorIDs/" + group + "-hugo.txt"
            elif species == "mouse":
                config_path_to_regs = __get_pyviper_dir() + "/data/regulatorIDs/" + group + "-hugo-mouse.txt"
        path_to_regs = config_path_to_regs
    regs_list = __load_regulators_from_txt(path_to_regs)
    return regs_list

# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTIONS ------------------------------
# -----------------------------------------------------------------------------

# Use can specify path: otherwise it goes to config default
def load_TFs(species = None, path_to_tfs = None):
    tfs_list = _load_regulators("tfs", species, path_to_tfs)
    return(tfs_list)

def load_coTFs(species = None, path_to_cotfs = None):
    cotfs_list = _load_regulators("cotfs", species, path_to_cotfs)
    return(cotfs_list)

def load_sig(species = None, path_to_sig = None):
    sig_list = _load_regulators("sig", species, path_to_sig)
    return(sig_list)

def load_surf(species = None, path_to_surf = None):
    surf_list = _load_regulators("surf", species, path_to_surf)
    return(surf_list)
