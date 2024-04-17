import pandas as pd
import numpy as np
import scanpy as sc
from ._filtering_funcs import *
from ._filtering_funcs import _get_anndata_filtered_by_feature_group
from ._viper import viper
from ._load._load import msigdb_regulon
from .interactome import Interactome
from ._pp import _nes_to_pval_df, _sig_clusters_adata

def _pca(adata,
         layer=None,
         filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
         **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)

    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]

def _dendrogram(adata,
               groupby,
               key_added=None,
               layer=None,
               filter_by_feature_groups=None, # ["tfs", "cotfs", "sig", "surf"],
               **kwargs):
    adata_filt = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups)
    if key_added is None:
        key_added = f'dendrogram_{groupby}'
    sc.tl.dendrogram(adata_filt, groupby, **kwargs, key_added = key_added)
    adata.uns[key_added] = adata_filt.uns[key_added]


def _generate_interactome_from_pax_data(pax_data,
                                        interactome_name="vpmat",
                                        n_top=50,
                                        is_symmetric=True):

    if isinstance(pax_data, anndata.AnnData):
        vpmat = pax_data.to_df()
    elif isinstance(pax_data, pd.DataFrame):
        vpmat = pax_data
    else:
        raise ValueError("pax_data must be anndata.AnnData or pd.DataFrame.")


    n_cells = vpmat.shape[0]
    n_mrs = vpmat.shape[1]
    cell_names = vpmat.index

    # For each sample, we calculate index arragement that would sort the vector
    sorted_order_array = np.argsort(-vpmat.values)
    # We then get the MRs ranked for each sample by indexing with this sorted order
    mrs_ranked_array = np.array(vpmat.columns)[sorted_order_array]

    if is_symmetric:
        # Slice the top n_top/2 rows and bottom n_top/2 rows
        # Get the top 25 and bottom 25 rows
        n_top_half = int(n_top/2)
        selected_column_indices = list(range(0,n_top_half)) + list(range(n_mrs-n_top_half,n_mrs))
        cell_i_mor = np.concatenate((np.ones(n_top_half), np.full(n_top_half, -1)))
    else:
        # Slice the top n_top rows
        selected_column_indices = list(range(0,n_top))
        cell_i_mor = np.ones(n_top)

    top_mrs_ranked_array = mrs_ranked_array[:, selected_column_indices]
    top_mrs_ranked_array_1D = top_mrs_ranked_array.flatten()

    regulator = np.repeat(cell_names, n_top)
    target = top_mrs_ranked_array_1D
    mor = np.tile(cell_i_mor, n_cells)

    net_table = pd.DataFrame({
            'regulator': regulator,
            'target': top_mrs_ranked_array_1D,
            'mor': mor,
            'likelihood': 1
        })

    return Interactome(interactome_name, net_table)

def _oncomatch(pax_data_to_test,
                pax_data_for_cMRs,
                tcm_size = 50,
                both_ways = False,
                om_max_NES_threshold = 30,
                om_min_logp_threshold = 0,
                lower_tail = True,
                enrichment = 'aREA',
                key_added = 'om',
                return_as_df = False,
                copy = False):
    if copy: pax_data_to_test = pax_data_to_test.copy()

    if enrichment is None:
        enrichment = 'narnea'
    else:
        enrichment = enrichment.lower()

    if isinstance(pax_data_to_test, anndata.AnnData):
        vpmat_to_test = pax_data_to_test.to_df()
    elif isinstance(pax_data_to_test, pd.DataFrame):
        vpmat_to_test = pax_data_to_test
    else:
        raise ValueError("pax_data_to_test must be anndata.AnnData or pd.DataFrame.")

    if isinstance(pax_data_for_cMRs, anndata.AnnData):
        vpmat_for_cMRs = pax_data_for_cMRs.to_df()
    elif isinstance(pax_data_for_cMRs, pd.DataFrame):
        vpmat_for_cMRs = pax_data_for_cMRs
    else:
        raise ValueError("pax_data_for_cMRs must be anndata.AnnData or pd.DataFrame.")

    # Compute intersection of regulons
    regs_in_common = np.intersect1d(vpmat_to_test.columns.values, vpmat_for_cMRs.columns.values)
    vpmat_to_test = vpmat_to_test[regs_in_common]
    vpmat_for_cMRs = vpmat_for_cMRs[regs_in_common]


    if both_ways is True:
        interactome_test = _generate_interactome_from_pax_data(
            vpmat_to_test,
            interactome_name = "vpmat_to_test",
            n_top = tcm_size
        )
        interactome_cMRs = _generate_interactome_from_pax_data(
            vpmat_for_cMRs,
            interactome_name = "vpmat_to_test",
            n_top = tcm_size
        )

        om_t = viper(gex_data=anndata.AnnData(vpmat_to_test, dtype='float64'),
                     interactome=interactome_cMRs,
                     enrichment=enrichment,
                     min_targets=0,
                     output_as_anndata=False,
                     verbose=False)
        if enrichment == 'narnea': om_t = om_t['nes']

        om_q = viper(gex_data=anndata.AnnData(vpmat_for_cMRs, dtype='float64'),
                     interactome=interactome_test,
                     enrichment=enrichment,
                     min_targets=0,
                     output_as_anndata=False,
                     verbose=False)
        if enrichment == 'narnea': om_q = om_q['nes']

        # Replace NaN (missing) values with 0 in om_t
        om_t[np.isnan(om_t)] = 0

        # Replace NaN (missing) values with 0 in om_q
        om_q[np.isnan(om_q)] = 0

        # Clip values greater than om_max_NES_threshold in om_t
        om_t = np.where(om_t > om_max_NES_threshold, om_max_NES_threshold, om_t)

        # Clip values greater than om_max_NES_threshold in om_q
        om_q = np.where(om_q > om_max_NES_threshold, om_max_NES_threshold, om_q)

        # Tranpose om_q so it has the same shape as om_t
        om_q = np.transpose(om_q)

        # Average them
        om = (om_t+om_q)/2
    else:
        interactome_cMRs = _generate_interactome_from_pax_data(
            vpmat_for_cMRs,
            interactome_name = "vpmat_to_test",
            n_top = tcm_size
        )
        om = viper(gex_data=anndata.AnnData(vpmat_to_test, dtype='float64'),
                   interactome=interactome_cMRs,
                   enrichment=enrichment,
                   min_targets=0,
                   output_as_anndata=False,
                   verbose=False)
        if enrichment == 'narnea': om = om['nes']

        # Replace NaN (missing) values with 0 in om
        om[np.isnan(om)] = 0

        # Clip values greater than om_max_NES_threshold in om
        om = np.where(om > om_max_NES_threshold, om_max_NES_threshold, om)

    # Compute p_values
    # cond = om < 7
    # om[cond] = 1 - norm.cdf(om[cond]) # accurate (i.e. same as R) when NES scores are small (returns 0.0 when given big NES, e.g. 10)
    # om[~cond] = norm.logcdf(om[~cond])*-1 # accurate (i.e. same as R) when NES scores are big (e.g. 10)

    om = pd.DataFrame(om, index = vpmat_to_test.index, columns = vpmat_for_cMRs.index)
    om = _nes_to_pval_df(om, lower_tail=lower_tail)

    # Log transform
    om = -np.log10(om)

    # Clip values smaller than om_min_logp_threshold in om
    om = np.where(om < om_min_logp_threshold, 0, om)

    om = pd.DataFrame(om, index = vpmat_to_test.index, columns = vpmat_for_cMRs.index)

    if return_as_df: return om

    pax_data_to_test.obsm[key_added] = om

    if copy: return pax_data_to_test

def _find_top_mrs_from_sig(stouffer_sig, N, both, rank):
    if rank is False:
        top_mrs_clusts_df = pd.DataFrame(False,index=stouffer_sig.columns, columns=stouffer_sig.index)
    else:
        top_mrs_clusts_df = pd.DataFrame(0,index=stouffer_sig.columns, columns=stouffer_sig.index)

    for i in range(stouffer_sig.shape[0]):
        stouffer_sig_clust_i = stouffer_sig.iloc[i,:]
        sorted_mrs_clust_i = stouffer_sig_clust_i.index.values[np.flip(np.argsort(stouffer_sig_clust_i.values))].flatten()
        top_mrs_clust_i = stouffer_sig_clust_i[sorted_mrs_clust_i[0:N]].index.values
        if both is True:
            bottom_mrs_clust_i = stouffer_sig_clust_i[sorted_mrs_clust_i[-N:]].index.values
            top_mrs_clust_i = np.concatenate((top_mrs_clust_i, bottom_mrs_clust_i))
        if rank is False:
            top_mrs_clusts_df.iloc[:, i] = np.isin(stouffer_sig_clust_i.index, top_mrs_clust_i)
        else:
            N_to_1 = np.flip(np.arange(N)+1)
            neg_1_to_neg_N = -1*(np.arange(N)+1)
            N_to_neg_N = np.concatenate((N_to_1, neg_1_to_neg_N))
            if both is True:
                top_mrs_clusts_df.iloc[:, i].loc[top_mrs_clust_i] = N_to_neg_N
            else:
                top_mrs_clusts_df.iloc[:, i].loc[top_mrs_clust_i] = N_to_1
    return top_mrs_clusts_df

def _find_top_mrs(adata,
                  pca_slot = "X_pca",
                  obs_column_name = None,
                  layer = None,
                  N = 50,
                  both = True,
                  method = "stouffer",
                  key_added = "mr",
                  filter_by_feature_groups=None,
                  rank=False,
                  filter_by_top_mrs = False,
                  return_as_df = False,
                  copy = False):
    if copy and return_as_df:
        raise ValueError("copy and return_as_df cannot both be True.")
    if copy: adata = adata.copy()

    sig = _sig_clusters_adata(adata,
                              obs_column_name,
                              layer,
                              filter_by_feature_groups,
                              sig_method = method,
                              compute_pvals = False,
                              pca_slot = "X_pca")
    result_df = _find_top_mrs_from_sig(sig, N, both, rank)
    result_df.columns.str.replace('_scores', '')

    if obs_column_name is None:
        result_df.columns = [key_added]
    else:
        result_df.columns = key_added + "_" + result_df.columns
    adata.var = pd.concat([adata.var, result_df], axis=1, join='inner')

    if filter_by_top_mrs:
        adata._inplace_subset_var(adata.var[key_added].values.flatten()==True)

    if return_as_df:
        mrs_df = result_df.apply(lambda col: result_df.index[col].tolist())
        return mrs_df
    elif copy:
        return adata

def _path_enr(gex_data,
             pathway_interactome,
             layer=None,
             eset_filter=True,
             method=None,  # [None, "scale", "rank", "mad", "ttest"],
             enrichment='aREA',  # [None, 'area','narnea'],
             mvws=1,
             njobs=1,
             batch_size=10000,
             verbose=True,
             output_as_anndata=True,
             transfer_obs=True,
             store_input_data=True
             ):
    if isinstance(pathway_interactome, str):
        collection = pathway_interactome.lower()
        if collection in ["c2", "c5", "c6", "c7", "h"]:
            pathway_interactome = msigdb_regulon(collection)
        else:
            raise ValueError(
                'pathway_interactome "' + str(pathway_interactome) + '" is not in "c2", "c5", "c6", "c7", "h".'
            )

    pathway_interactome.filter_targets(gex_data.var_names)
    return viper(
        gex_data,
        pathway_interactome,
        layer,
        eset_filter,
        method,
        enrichment,
        mvws,
        0, #min_targets
        njobs,
        batch_size,
        verbose,
        output_as_anndata,
        transfer_obs,
        store_input_data
    )
