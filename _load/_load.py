### Import dependencies
from ._load_translate import *
from ._load_regulators import *
from ._load_msigdb import *

### EXPORT LIST
__all__ = ['mouse2human', 'human2mouse',
           'TFs', 'coTFs', 'sig', 'surf',
           'msigdb_regulon']

def mouse2human():
    return load_mouse2human()

def human2mouse():
    return load_human2mouse()

def TFs(species = None, path_to_tfs = None):
    return load_TFs(species, path_to_tfs)

def coTFs(species = None, path_to_cotfs = None):
    return load_coTFs(species, path_to_cotfs)

def sig(species = None, path_to_sig = None):
    return load_sig(species, path_to_sig)

def surf(species = None, path_to_surf = None):
    return load_surf(species, path_to_surf)

def msigdb_regulon(collection = "h"):
    return load_msigdb_regulon(collection)
