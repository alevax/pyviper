import numpy as np
from scipy.stats import norm
from scipy.stats import rankdata
import pandas as pd
import anndata
import scanpy as sc
from pyther_classes import *
import pathlib

def get_pyther_dir():
    pyther_dir = str(pathlib.Path(__file__).parent.parent.resolve())
    return(pyther_dir)

def interactome_from_tsv(filePath, intName):
    """\
    Allows the user to load an interactome object from a TSV file.

    The TSV file is created by the R function InteractomeToTable.

    Parameters
    ----------
    filePath
        The file to the regulon.tsv file.
    intName
        The name of the interactome.
    Returns
    -------
    A dictionary of :class:`~pyther_classes.Interactome`.
    """
    # read file
    netTable = pd.read_csv(filePath, sep = '\t')
    intObj = Interactome('intName')
    # loop through regulators
    uniqueRegs = netTable.regulator.unique()
    for u in uniqueRegs:
        # subset dataframe
        uDF = netTable[netTable.regulator == u]
        # make dictionaries
        icDict = dict(zip(uDF.target, uDF.likelihood))
        morDict = dict(zip(uDF.target, uDF.mor))
        # make regulon object
        regObj = Regulon(u, icDict, morDict)
        intObj.addReg(u, regObj)
    # return
    return(intObj)

def aREA(gex_data, intObj, layer = None):
    """\
    Allows the individual to infer normalized enrichment scores from gene
    expression data using the analytical ranked enrichment analysis (aREA)
    function.

    It is the basis of the VIPER (Virtual Inference of Protein-activity
    by Enriched Regulon analysis) algorithm.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object (e.g. from Scanpy).
    intObj
        The interactome object.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
    """
    if layer is None:
        gesMat = gex_data.X
    else:
        gesMat = gex_data.layers[layer]

    # rank transform the GES
    rankMat = rankdata(gesMat, axis = 1)

    # find intersecting genes
    targetSet = intObj.get_targetSet()
    varNames = gex_data.var_names.to_list()
    intersectGenes = [value for value in targetSet if value in varNames]

    # reduce regulon matrices
    icMat = intObj.icMat().loc[intersectGenes]
    morMat = intObj.morMat().loc[intersectGenes]

    # prepare the 1-tailed / 2-tailed matrices
    gesInds = [varNames.index(i) for i in intersectGenes]
    ges2T = rankMat / (rankMat.shape[1] + 1)
    ges1T = abs(ges2T - 0.5) * 2
    ges1T = ges1T + (1 - np.max(ges1T))/2
    ges2TQ = norm.ppf(ges2T[:, gesInds])
    ges1TQ = norm.ppf(ges1T[:, gesInds])

    # 2-tail enrichment
    dES = pd.DataFrame.transpose(pd.DataFrame.mul(icMat, morMat))
    dES = dES.dot(np.transpose(ges2TQ))

    # 1-tail enrichemnt
    uES = pd.DataFrame.transpose(pd.DataFrame.mul(1 - abs(morMat), icMat))
    uES = uES.dot(np.transpose(ges1TQ))

    # integrate
    iES = (abs(dES) + uES * (uES > 0)) * np.sign(dES)
    # make NES
    nES = iES.mul(intObj.icpVec(), 0)
    nES = np.transpose(nES)
    nES.index = gex_data.obs.index

    return(nES)

def pyther(gex_data, intObj, layer = None):
    """\
    Allows the individual to infer protein activity from gene expression using
    the VIPER (Virtual Inference of Protein-activity by Enriched Regulon
    analysis) algorithm.

    Parameters
    ----------
    gex_data
        Gene expression stored in an anndata object.
    intObj
        The interactome object.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~anndata._core.anndata.AnnData` with NES values
    stored in the .X slot and the original anndata object with gene expression
    stored in the .gex_data slot.
    """
    # aREA takes gex_data.X
    nesMat = aREA(gex_data, intObj, layer)
    # Create an Anndata object from the nesMat
    pax_data = mat_to_anndata(nesMat)
    # Store the GExpr Anndata object in the PAct Anndata object
    pax_data.gex_data = gex_data
    return(pax_data)

def path_enr(adata, intObj, layer = None):
    """\
    Allows the individual to infer normalized enrichment scores of pathways
    using the analytical ranked enrichment analysis (aREA) function.

    This is a wrapper for the aREA function in Pyther.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    intObj
        The interactome object containing pathways.
    layer
        The layer in the anndata object to use as the gene expression input
        (default = None).
    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame` containing NES values.
    """
    # aREA takes the pathways interactome and the adata
    pathEnrMat = aREA(adata, intObj, layer)
    # Create a new Anndata object
    pathEnrObj = mat_to_anndata(pathEnrMat)
    # This means we did pathway enrichment on VIPER: adata is pax_data
    if adata.gex_data is not None:
        pathEnrObj.gex_data = adata.gex_data
        adata.gex_data = None
        pathEnrObj.pax_data = adata
    # This means we did pathway enrichment on gExpr: adata is gex_data
    else:
        pathEnrObj.gex_data = adata
    return(pathEnrObj)

def mat_to_anndata(mat):
    # Helper function for *pyther* and *path_enr*
    # Create obs dataframe
    mat_sampleNames = pd.DataFrame(index=range(len(mat.index.values)),columns=range(0))
    mat_sampleNames.index = mat.index.values
    mat_sampleNames

    # Create var dataframe
    mat_features = pd.DataFrame(index=range(len(mat.columns.values)),columns=range(0))
    mat_features.index = mat.columns.values
    mat_features

    # Convert the pandas dataframe from Pyther into a new Anndata object
    pax_data = anndata.AnnData(X=mat,
                               obs=mat_sampleNames,
                               var=mat_features)
    return(pax_data)

def tl_pca(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, **kwargs)
    adata.obsm["X_pca"] = adata_filt.obsm["X_pca"]
    return(adata)

def tl_tsne(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.tsne(adata_filt, **kwargs)
    adata.obsm["X_tsne"] = adata_filt.obsm["X_tsne"]
    return(adata)

def tl_umap(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            svd_solver='arpack',
            n_neighbors=10,
            n_pcs=40,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    # Create a second anndata object
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # Move sc.pp.neighbors results to adata_filt
    # adata_filt.obsp["distances"] = adata.obsp["distances"]
    # adata_filt.obsp["connectivities"] = adata.obsp["connectivities"]
    # adata_filt.uns["neighbors"] = adata.uns["neighbors"]
    # Compute the UMAP
    sc.tl.umap(adata_filt, **kwargs)
    # Give the UMAP from the second anndata object to the original
    adata.obsm["X_umap"] = adata_filt.obsm["X_umap"]
    return(adata)

#KeyError: 'No "neighbors" in .uns'
def tl_draw_graph(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.draw_graph(adata_filt, **kwargs)
    if "X_draw_graph_fa" in list(pAct_adata.obsm.keys()):
        adata.obsm["X_draw_graph_fa"] = adata_filt.obsm["X_draw_graph_fa"]
    if "X_draw_graph_fr" in list(pAct_adata.obsm.keys()):
        adata.obsm["X_draw_graph_fr"] = adata_filt.obsm["X_draw_graph_fr"]
    return(adata)

#ValueError: You need to run `pp.neighbors` first to compute a neighborhood graph.
def tl_diffmap(adata,
            *,
            filter_by_feature_groups=None, #["TFs", "CoTFs", "sig", "surf"],
            path_to_tfs = None,
            path_to_cotfs = None,
            path_to_sig = None,
            path_to_surf = None,
            **kwargs):
    if filter_by_feature_groups is None:
        filter_by_feature_groups = "all"
    adata_filt = get_anndata_filtered_by_feature_group(adata,
                                                       filter_by_feature_groups,
                                                       path_to_tfs,
                                                       path_to_cotfs,
                                                       path_to_sig,
                                                       path_to_surf)
    sc.tl.pca(adata_filt, svd_solver=svd_solver)
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.diffmap(adata_filt, **kwargs)
    adata.obsm["X_diffmap"] = adata_filt.obsm["X_diffmap"]
    adata.uns['diffmap_evals'] = adata_filt.uns['diffmap_evals']
    return(adata)

def get_anndata_filtered_by_feature_group(adata,
                               feature_groups="all", #["TFs", "CoTFs", "sig", "surf"],
                               path_to_tfs = None,
                               path_to_cotfs = None,
                               path_to_sig = None,
                               path_to_surf = None,):
    features_list = get_features_list(adata,
            feature_groups,
            path_to_tfs,
            path_to_cotfs,
            path_to_sig,
            path_to_surf)
    adata_filt = get_anndata_filtered_by_feature_list(adata, features_list)
    return(adata_filt)

def get_anndata_filtered_by_feature_list(adata, features_list):
    features_indices = get_features_indices(adata, features_list)
    mat = get_mat_from_anndata(adata, features_indices)
    adata_with_features_only = mat_to_anndata(mat)
    return(adata_with_features_only)

def get_mat_from_anndata(adata, features_indices):
    mat = pd.DataFrame(adata.X[:,features_indices],
                       index = list(adata.obs_names),
                       columns = adata.var_names[features_indices,])
    return(mat)

def get_features_indices(adata, features_list):
    true_false_list = pd.Series(adata.var_names).isin(features_list).tolist()
    features_indices = [i for i, x in enumerate(true_false_list) if x]
    return(features_indices)

def get_features_list(adata,
                    feature_groups="all",
                    path_to_tfs = None,
                    path_to_cotfs = None,
                    path_to_sig = None,
                    path_to_surf = None):
    if(type(feature_groups) is str):
        feature_groups = [feature_groups]
    if "all" not in feature_groups:
        feature_groups = [x.lower() for x in feature_groups]
        features_list = list()
        if "tfs" in feature_groups or "tf" in feature_groups:
            tfs = load_TFs(path_to_tfs)
            features_list.extend(tfs)
        if "cotfs" in feature_groups or "cotf" in feature_groups:
            cotfs = load_coTFs(path_to_cotfs)
            features_list.extend(cotfs)
        if "sigs" in feature_groups or "sig" in feature_groups:
            sig = load_sig(path_to_sig)
            features_list.extend(sig)
        if "surfs" in feature_groups or "surf" in feature_groups:
            surf = load_surf(path_to_surf)
            features_list.extend(surf)
        features_list = intersection(features_list, list(adata.var_names))
    else:
        features_list = adata.var_names
    return(features_list)

def load_TFs(path_to_tfs = None):
    if path_to_tfs is None:
        path_to_tfs = get_pyther_dir() + "/data/regulatorIDs/tfs-hugo.txt"
    tfs_list = load_regulators(path_to_tfs)
    return(tfs_list)

def load_coTFs(path_to_cotfs = None):
    if path_to_cotfs is None:
        path_to_cotfs = get_pyther_dir() + "/data/regulatorIDs/cotfs-hugo.txt"
    cotfs_list = load_regulators(path_to_cotfs)
    return(cotfs_list)

def load_sig(path_to_sig = None):
    if path_to_sig is None:
        path_to_sig = get_pyther_dir() + "/data/regulatorIDs/sig-hugo.txt"
    sig_list = load_regulators(path_to_sig)
    return(sig_list)

def load_surf(path_to_surf = None):
    if path_to_surf is None:
        path_to_surf = get_pyther_dir() + "/data/regulatorIDs/surface-hugo.txt"
    surf_list = load_regulators(path_to_surf)
    return(surf_list)

def load_regulators(path_to_txt):
    with open(path_to_txt) as temp_file:
        regulator_set = [line.rstrip('\n') for line in temp_file]
    return(regulator_set)

def load_msigdb_regulon(collection = "c2"):
    reg = None
    if(collection.lower() == "c2"):
        reg_path = get_pyther_dir() + "/data/regulons/msigdb-c2-as-regulon.tsv"
        reg = interactome_from_tsv(reg_path, "MSigDB_C2")
    elif(collection.lower() == "c5"):
        reg_path = get_pyther_dir() + "/data/regulons/msigdb-c5-as-regulon.tsv"
        reg = interactome_from_tsv(reg_path, "MSigDB_C5")
    elif(collection.lower() == "c6"):
        reg_path = get_pyther_dir() + "/data/regulons/msigdb-c6-as-regulon.tsv"
        reg = interactome_from_tsv(reg_path, "MSigDB_C6")
    elif(collection.lower() == "c7"):
        reg_path = get_pyther_dir() + "/data/regulons/msigdb-c7-as-regulon.tsv"
        reg = interactome_from_tsv(reg_path, "MSigDB_C7")
    elif(collection.lower() == "h"):
        reg_path = get_pyther_dir() + "/data/regulons/msigdb-h-as-regulon.tsv"
        reg = interactome_from_tsv(reg_path, "MSigDB_H")
    return(reg)

# Python program to illustrate the intersection
# of two lists in most simple way
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def pl_umap(adata,
            *,
            plot_stored_gex_data=False,
            plot_stored_pax_data=False,
            **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.umap.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.umap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.umap(adata, **kwargs)

def pl_scatter(adata,
               *,
               plot_stored_gex_data=False,
               plot_stored_pax_data=False,
               **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.scatter.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.scatter function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.scatter(adata,**kwargs)

def pl_heatmap(adata,
     *,
     plot_stored_gex_data=False,
     plot_stored_pax_data=False,
     **kwargs):
    """\
    A wrapper for the scanpy function sc.pl.heatmap.

    Parameters
    ----------
    adata
        Gene expression or protein activity stored in an anndata object.
    plot_stored_gex_data
        Plot stored gex_data on the protein activity (or pathways) UMAP.
    plot_stored_pax_data
        If the adata is the output of the path_enr function and pax_data is
        stored in adata, this allows one to plot pax_data values in the UMAP
        created using pathways results.
    **kwargs
        Arguments to provide to the sc.pl.heatmap function.
    Returns
    -------
    A plot of :class:`~matplotlib.axes.Axes`.
    """
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.heatmap(adata,**kwargs)

def pl_dotplot(adata,
                      *,
                      plot_stored_gex_data = False,
                      plot_stored_pax_data = False,
                      **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.dotplot(adata,**kwargs)

def pl_tracksplot(adata,
                      *,
                      plot_stored_gex_data = False,
                      plot_stored_pax_data = False,
                      **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.tracksplot(adata,**kwargs)

def pl_violin(adata,
                     *,
                     plot_stored_gex_data = False,
                     plot_stored_pax_data = False,
                     **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.violin(adata,**kwargs)

def pl_stacked_violin(adata,
                             *,
                             plot_stored_gex_data = False,
                             plot_stored_pax_data = False,
                             **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.stacked_violin(adata,**kwargs)

def pl_matrixplot(adata,
                         *,
                         plot_stored_gex_data = False,
                         plot_stored_pax_data = False,
                         **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.matrixplot(adata,**kwargs)

def pl_clustermap(adata,
                         *,
                         plot_stored_gex_data = False,
                         plot_stored_pax_data = False,
                         **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.clustermap(adata,**kwargs)

def pl_ranking(adata,
                      *,
                      plot_stored_gex_data = False,
                      plot_stored_pax_data = False,
                      **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.ranking(adata,**kwargs)

def pl_dendrogram(adata,
                         *,
                         plot_stored_gex_data = False,
                         plot_stored_pax_data = False,
                         **kwargs):
    if(plot_stored_gex_data is True):
        adata = get_gExpr_anndata_with_NES_umap(adata)
    elif(plot_stored_pax_data is True):
        adata = get_pAct_anndata_with_pathEnr_umap(adata)
    sc.pl.dendrogram(adata,**kwargs)

def get_gExpr_anndata_with_NES_umap(adata):
    adata_gExpr = adata.gex_data
    adata_gExpr.obsm["X_umap"] = adata.obsm["X_umap"]
    return(adata_gExpr)

def get_pAct_anndata_with_pathEnr_umap(adata):
    adata_pAct = adata.pax_data
    adata_pAct.obsm["X_umap"] = adata.obsm["X_umap"]
    return(adata_pAct)
