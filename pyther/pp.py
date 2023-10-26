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


def aracne3_to_regulon(net_file, net_df=None, anno=None, MI_thres=0, regul_size=50, is_for_viper=True, with_count_values=False,
                       to_dataframe = True):
    if net_df is None:
        net = pd.read_csv(net_file, sep='\t')
    else:
        net = net_df.copy()

    if anno is not None:
        if anno.shape[1] != 2 or not isinstance(anno, (pd.DataFrame, pd.Series, pd.Matrix)):
            raise ValueError("anno should contain two columns: 1-original symbol, 2-new symbol")

        anno.columns = ["old", "new"]

        ## Convert gene symbols:
        net['regulator.values'] = anno.set_index('old').loc[net['regulator.values'], 'new'].values
        net['target.values'] = anno.set_index('old').loc[net['target.values'], 'new'].values

    ## Network filtering
    net = net[net['mi.values'] > MI_thres]

    ## Total MR set
    mr = net['regulator.values'].unique()
    if to_dataframe:
        regul = pd.DataFrame(columns=["regulator","target","mor","likelihood"])
        for mri in mr:
            tmp_net_data = net[net['regulator.values'] == mri].copy()
            ## Sort interactions in decreasing order of count and MI
            tmp_net_data.sort_values(by=['count.values', 'mi.values'], ascending=[False, False], inplace=True)

            # print(tmp_net_data.shape)
            ## Top 50 interactions
            tmp_net_data = tmp_net_data.iloc[:min(regul_size, len(tmp_net_data))]

            ## Regulatory mode = spearman correlation score
            tmp_net_data['am.values'] = tmp_net_data['scc.values']

            if with_count_values:
                tmp_net_data['count.values'] = tmp_net_data['count.values']

            ## Regulatory weight = scaled MI
            tmp_net_data['aw.values'] = tmp_net_data['mi.values'] / tmp_net_data['mi.values'].max()

            for index, row in tmp_net_data.iterrows():
                regul.loc[len(regul.index)] = row.loc[["regulator.values", "target.values","aw.values","am.values"]].array
    else:
        regul = {}
        for mri in mr:
            ## Regulon for MR_i
            tmp_net_data = net[net['regulator.values'] == mri].copy()
            ## Sort interactions in decreasing order of count and MI
            tmp_net_data.sort_values(by=['count.values', 'mi.values'], ascending=[False, False], inplace=True)

            # print(tmp_net_data.shape)
            ## Top 50 interactions
            tmp_net_data = tmp_net_data.iloc[:min(regul_size, len(tmp_net_data))]

            ## Regulatory mode = spearman correlation score
            tmp_net_data['am.values'] = tmp_net_data['scc.values']

            if with_count_values:
                tmp_net_data['count.values'] = tmp_net_data['count.values']

            ## Regulatory weight = scaled MI
            tmp_net_data['aw.values'] = tmp_net_data['mi.values'] / tmp_net_data['mi.values'].max()
            if is_for_viper:
                if with_count_values:
                    tmp_regul = {'likelihood': {}, 'tfmode': {}, 'subnets':{}}
                    for index, row in tmp_net_data.iterrows():
                    
                        tmp_regul['likelihood'][row["target.values"]] = row['aw.values']
                        tmp_regul['tfmode'][row["target.values"]] = row['am.values']
                        tmp_regul['subnets'][row["target.values"]] = row['count.values']
                        
                else:
                    tmp_regul = {'likelihood': {}, 'tfmode': {}}
                    for index, row in tmp_net_data.iterrows():
                        tmp_regul['likelihood'][row["target.values"]] = row['aw.values']
                        tmp_regul['tfmode'][row["target.values"]] = row['am.values']
                    
                    
            else:
                tmp_regul = {'am': {}, 'aw': {}}
                for index, row in tmp_net_data.iterrows():
                    tmp_regul['am'][row["target.values"]] = row['aw.values']
                    tmp_regul['aw'][row["target.values"]] = row['am.values']

                    

            regul[mri] = tmp_regul
    return regul