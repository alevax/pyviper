
### Import dependencies
from ._load_get_path import __get_pyviper_dir
import pandas as pd

### EXPORT LIST
# __all__ = ['load_mouse2human', 'load_human2mouse']
__all__ = ['load_human2mouse']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------ LOAD TRANSLATE DFS FUNCTIONS -----------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# def load_mouse2human():
#     mouse2human = pd.read_csv(__get_pyviper_dir() + "/data/translate/human2mouse.csv", dtype=str)
#     del mouse2human[mouse2human.columns[0]]
#     return(mouse2human)
# def load_human2mouse():
#     human2mouse = pd.read_csv(__get_pyviper_dir() + "/data/translate/mouse2human.csv", dtype=str)
#     del human2mouse[human2mouse.columns[0]]
#     return(human2mouse)

def load_human2mouse():
    human2mouse = pd.read_csv(__get_pyviper_dir() + "/data/translate/human_mouse_gene_translation.csv", dtype=str)
    # del human2mouse[human2mouse.columns[0]]
    return(human2mouse)
