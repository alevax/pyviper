
### Import dependencies
from ._load_get_path import __get_vithon_dir
import pandas as pd

### EXPORT LIST
__all__ = ['load_mouse2human', 'load_human2mouse']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------ LOAD TRANSLATE DFS FUNCTIONS -----------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def load_mouse2human():
    mouse2human = pd.read_csv(__get_vithon_dir() + "/data/translate/human2mouse.csv")
    del mouse2human[mouse2human.columns[0]]
    return(mouse2human)
def load_human2mouse():
    human2mouse = pd.read_csv(__get_vithon_dir() + "/data/translate/mouse2human.csv")
    del human2mouse[human2mouse.columns[0]]
    return(human2mouse)
