### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import anndata
from scipy.stats import rankdata

### ---------- EXPORT LIST ----------
__all__ = []


def _mad_from_R(x, center=None, constant=1.4826, low=False, high=False):
    if center is None:
        center=np.median(x)
    x = x[~np.isnan(x)] if np.isnan(x).any() else x
    n = len(x)
    if (low or high) and n % 2 == 0:
        if low and high:
            raise ValueError("'low' and 'high' cannot be both True")
        n2 = n // 2 + int(high)
        return constant * np.sort(np.abs(x - center))[n2]
    return constant * np.median(np.abs(x - center))

# Function assumes features as rows and observations as columns
# Numerator Functions:
    # Median - numpy.median
    # Mean - numpy.mean
# Denominator Functions:
    # Median absolute deviation - mad_from_R
    # Standard deviation - statistics.stdev
def _rank_norm(x, NUM_FUN=np.median, DEM_FUN = _mad_from_R, verbose = False):
    rank = rankdata(x, axis=0)
    median = NUM_FUN(rank, axis=1, keepdims=True)#np.median(rank, axis=1, keepdims=True)
    mad = np.apply_along_axis(DEM_FUN, 1, rank)

    x = ((rank - median)/mad[:, np.newaxis])

    if verbose: print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
    if verbose: print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
    if verbose: print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))
    if verbose: print("- Features to Remove:")

    # Take care of infinite values
    max_finite = np.nanmax(x[np.isfinite(x)])
    min_finite = np.nanmin(x[np.isfinite(x)])
    x[np.isposinf(x)] = max_finite
    x[np.isneginf(x)] = min_finite

    x = np.where(np.isnan(x), np.nanmin(x), x)
    x = np.clip(x, a_min=np.nanmin(x), a_max=np.nanmax(x))
    if verbose: print("- Removing NULL/NA features ...")
    x = x[~np.isnan(x).any(axis=1)]

    if verbose: print("- Number of NA features:", np.sum(np.sum(np.isnan(x), axis=1)))
    if verbose: print("- Number of Inf features:", np.sum(np.isinf(np.sum(x, axis=1))))
    if verbose: print("- Number of 0 features:", np.sum(np.sum(x, axis=1) == 0))

    return x

def rank_norm(x, NUM_FUN=np.median, DEM_FUN = _mad_from_R, layer_input = None, layer_output = None, verbose = False):
    """\
    Compute a double rank normalization on an anndata, np.array, or pd.DataFrame.

    Parameters
    ----------
    x
        Data stored in an anndata object, np.array or pd.DataFrame.
    NUM_FUN (default: np.median)
        The first function to be applied across each column.
    DEM_FUN (default: _mad_from_R)
        The second function to be applied across each column.
    layer_input (default: None)
        For an anndata input, the layer to use. When None, the input layer is anndata.X.
    layer_output (default: None)
        For an anndata input, the name of the layer where to store.
        When None, this is anndata.X.
    verbose (default: False)
        Whether to give additional output about the progress of the function.

    Returns
    -------
    A double rank transformed version of the input data
    """
    if(isinstance(x, anndata.AnnData) or isinstance(x, anndata._core.anndata.AnnData)):
        if(layer_input is None):
            gesMat = x.X.copy()
        else:
            if(verbose): print('- Using the layer "' + layer_input + '" as input...')
            gesMat = x.layers[layer_input].copy()
    elif(isinstance(x, np.ndarray)):
        gesMat = x.copy()
    elif(isinstance(x, pd.DataFrame)):
        gesMat = x.copy().to_numpy
    else:
        raise Exception("In RankNorm(x), x must be anndata.AnnData, numpy.ndarray or pandas.DataFrame.")
    gesMat = _rank_norm(gesMat, NUM_FUN, DEM_FUN, verbose)
    if(isinstance(x, anndata.AnnData)):
        if(layer_output is None):
            x.X = gesMat
        else:
            if(verbose): print('- Saving result in the layer "' + layer_output + '"...')
            x.layers[layer_output] = gesMat
    else:
        x = gesMat
    return(x)
