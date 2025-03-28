### Import dependencies
import pathlib

### EXPORT LIST
__all__ = ['__get_pyviper_dir']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ GENERAL FUNCTION -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def __get_pyviper_dir():
    pyviper_dir = str(pathlib.Path(__file__).parent.parent.resolve())
    return(pyviper_dir)
