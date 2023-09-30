### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
from ._helpers import *
from .aREA.aREA_meta import aREA
from .NaRnEA.NaRnEA_meta import NaRnEA

### ---------- EXPORT LIST ----------
__all__ = ['pyther']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def sample_ttest(i,array):
    return ttest_1samp((array[i] - np.delete(array, i, 0)), 0).statistic

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def pyther(gex_data,
           interactome,
           layer=None,
           njobs=1,
           eset_filter=True,
           # bootstrap = 0,
           # dnull = None,
           # pleiotropy = False,
           # minsize=25,
           # adaptive_size=False,
           method=None,  # [None, "scale", "rank", "mad", "ttest"],
           enrichment='area',  # [None, 'area','narnea'],
           mvws=1,
           # pleiotropyArgs={'regulator':0.05, 'shadow':0.05, 'targets':10, "penalty":20, "method":"adaptive"},
           verbose=True,
           output_type='anndata',  # ['anndata', 'ndarray'],
           transfer_obs=True,
           store_gex_data=True
           ):
    gex_data_original = gex_data
    gex_data = gex_data_original.copy()

    pd.options.mode.chained_assignment = None

    if verbose:
        print("Preparing the association scores")

    ''' not in current version 3.8 of python , only in python 3.10
    match method:

        case 'scale': #scale each gene in all sample
            gex_data.X = (gex_data.X - np.mean(gex_data.X,axis=0))/np.std(gex_data.X,axis=0)
        case 'rank':
            gex_data.X = rankdata(gex_data.X,axis=0)*(np.random.random(gex_data.X.shape)*2/10-0.1)
        case 'mad':
            median = np.median(gex_data.X,axis=0)
            gex_data.X = (gex_data.X-median)/np.median(np.abs(gex_data.X-median))
        case 'ttest':
            # I dont quite understand this part maybe ask later:
            gex_data.X = gex_data.X
    '''
    if method == 'scale':
        gex_data.X = (gex_data.X - np.mean(gex_data.X,axis=0))/np.std(gex_data.X,axis=0)
    elif method == 'rank':
        gex_data.X = rankdata(gex_data.X,axis=0)*(np.random.random(gex_data.X.shape)*2/10-0.1)
    elif method == 'mad':
        median = np.median(gex_data.X,axis=0)
        gex_data.X = (gex_data.X-median)/np.median(np.abs(gex_data.X-median))
    elif method == 'ttest':
        gex_data.X = np.array([sample_ttest(i, gex_data.X.copy()) for i in range(gex_data.shape[0])])

#     if bootstrap > 0:
#
# #        num_sample = gex_data.X.shape[0]
#         bmean = np.zeros((bootstrap,gex_data.X.shape[1]))
#         bsd = np.zeros((bootstrap,gex_data.X.shape[1]))
#
#
#         sampled_indices = np.random.choice(gex_data.obs.index,bootstrap*len(gex_data.obs))
#
#         for i in range(bootstrap):
#             sample_name = sampled_indices[i*len(gex_data.obs):(i+1)*len(gex_data.obs)]
#             #sample mean
#             bmean[i] = np.mean(gex_data[sample_name,:].X, axis=0)
#             bsd[i] = np.std(gex_data[sample_name,:].X, axis=0)
#
#         #targets = intObj.get_targetSet()
#         # may need a bootstrap area
#         netMets = Parallel(n_jobs = njobs)(
#         (delayed)(bootstrap_aREA)(gex_data,iObj,eset_filter = True,layer = layer)
#         for iObj in interactome
#         )
#
#         preOp = consolidate_meta_aREA_results(netMets, mvws, verbose)
#
#     else:

    if enrichment is None: enrichment = 'narnea'
    if enrichment == 'area':
        if verbose: print("Computing regulons enrichment with aREA")
        preOp = aREA(gex_data, interactome, eset_filter, layer=layer,
                     mvws=mvws, njobs=njobs, verbose=verbose)
    elif enrichment == 'narnea':
        if verbose: print("Computing regulons enrichment with NaRnEa")
        preOp = NaRnEA(gex_data, interactome, layer=layer,
                       sample_weight=True, njobs=njobs, verbose=verbose)

    else:
        raise ValueError("Unsupported enrichment type:" + str(enrichment))

    if output_type == 'ndarray':
        op = preOp
    elif output_type == 'anndata':

        if enrichment == 'area':
            op = mat_to_anndata(preOp)
        else: #enrichment == 'narnea':
            op = mat_to_anndata(preOp["nes"])
            op.layers['pes'] = preOp["pes"]

        if transfer_obs is True:
            op.obs = op.obs.join(gex_data.obs)
        if store_gex_data is True:
            op.gex_data = gex_data_original
    else:
        raise ValueError("Unsupported output_type:" + str(output_type))

    return op #final result
