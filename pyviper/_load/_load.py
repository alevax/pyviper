### Import dependencies
from ._load_translate import *
from ._load_regulators import *
from ._load_msigdb import *
from ._load_pisces import *

### EXPORT LIST
__all__ = [#'mouse2human',
           'human2mouse',
           'TFs', 'coTFs', 'sig', 'surf',
           'msigdb_regulon',
           'pisces_network']

def human2mouse():
    """\
    Retrieves the human to mouse translation pd.DataFrame from pyVIPER's data
    folder. This dataframe contains six columns:
    human_symbol, mouse_symbol, human_ensembl, mouse_ensembl, human_entrez, mouse_entrez

    Returns
    -------
    A dataframe of :class:`~pandas.core.frame.DataFrame`.
    """
    return load_human2mouse()

def TFs(species = None, path_to_tfs = None):
    """\
    Retrieves a list of transcription factors (TFs).

    Parameters
    ----------
    species : default: None
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_tfs : default: None
        When left as None, the path to TFs setting in pyviper.config will be used.
        Otherwise, manually specify a filepath to a .txt file containing TFs,
        one on each line.

    Returns
    -------
    A list containing transcription factors.
    """
    return load_TFs(species, path_to_tfs)

def coTFs(species = None, path_to_cotfs = None):
    """\
    Retrieves a list of co-transcription factors (coTFs).

    Parameters
    ----------
    species : default: None
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_cotfs : default: None
        When left as None, the path to coTFs setting in pyviper.config will be used.
        Otherwise, manually specify a filepath to a .txt file containing coTFs,
        one on each line.

    Returns
    -------
    A list containing co-transcription factors.
    """
    return load_coTFs(species, path_to_cotfs)

def sig(species = None, path_to_sig = None):
    """\
    Retrieves a list of signalling proteins (sig).

    Parameters
    ----------
    species : default: None
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_sig : default: None
        When left as None, the path to sig setting in pyviper.config will be used.
        Otherwise, manually specify a filepath to a .txt file containing
        signaling proteins, one on each line.

    Returns
    -------
    A list containing signaling proteins.
    """
    return load_sig(species, path_to_sig)

def surf(species = None, path_to_surf = None):
    """\
    Retrieves a list of surface proteins (surf).

    Parameters
    ----------
    species : default: None
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_sig : default: None
        When left as None, the path to surf setting in pyviper.config will be used.
        Otherwise, manually specify a filepath to a .txt file containing surface
        proteins, one on each line.

    Returns
    -------
    A list containing signaling proteins.
    """
    return load_surf(species, path_to_surf)

def msigdb_regulon(collection):
    """
    Retrieves an object or a list of objects of class Interactome from pyviper's
    data folder containing a set of pathways from the Molecular Signatures
    Database (MSigDB), downloaded from https://www.gsea-msigdb.org/gsea/msigdb.
    
    Collections can be one of the following:

    - 'h' for Hallmark gene sets. Coherently expressed signatures derived by
       aggregating many MSigDB gene sets to represent well-defined biological
       states or processes.
    - 'c2' for curated gene sets. From online pathway databases, publications
       in PubMed, and knowledge of domain experts.
    - 'c5' for ontology gene sets. Consists of genes annotated by the same
       ontology term.
    - 'c6' for oncogenic signature gene sets. Defined directly from microarray
       gene expression data from cancer gene perturbations.
    - 'c7' for immunologic signature gene sets. Represents cell states and
        perturbations within the immune system.

    Parameters
    ----------
    collection : str or list of str
        A individual string or a list of strings containing the following:
        ["h", "c2", "c5", "c6", "c7"], corresponding to the collections above.

    Returns
    -------
    An individual object or list of objects of class pyviper.interactome.Interactome.
    """
    return load_msigdb_regulon(collection)

def pisces_network(tissue, celltype, desired_format = "human_ensembl"):
    """\
    The pipeline for Protein Activity Inference in Single Cells (PISCES) is a
    regulatory-network-based methdology for the analysis of single cell gene
    expression profiles. PISCES leverages the assembly of lineage-specific gene
    regulatory networks, to accurately measure activity of each protein based
    on the expression its transcriptional targets (regulon), using the ARACNe
    and metaVIPER algorithms, respectively. These networks are available to
    be loaded from the GitHub repository with this function. Reference using
    the citation below.

    Parameters
    ----------
    tissue
        The tissue type of the network. See options below.
    celltype
        The celltype of the network. See options below.
    desired_format : default: human_ensembl
        The gene name format of the regulators and targets in the network.
        Choose from: human_ensembl (default), human_symbol, human_entrez,
        mouse_ensembl, mouse_symbol, or mouse_entrez.

    Returns
    -------
    An individual object of class pyviper.interactome.Interactome.

    Overview of celltypes
    -------
    Choose from one of the following networks (tissue: celltype):
        adipose_tissue:
            adipocytes
            dendritic
            fibroblasts
            lymphoid
            myeloid
            smc
        bone_marrow
            lymphoid
        brain
            neuron
        breast
            adipocytes
            dendritic
            endothelial
            epithelial
            fibroblasts
            lymphoid
            myeloid
            smc
        bronchus
            epithelial
            keratinocytes
        colon
            enterocyte
            goblet
            lymphoid
        endometrium
            adipocytes
            endothelial
            epithelial
            fibroblasts
            lymphoid
        esophagus
            epithelial
            fibroblasts
        eye
            bipolar
            muller-glia
            photoreceptor
        heart_muscle
            cardiomyocyte
            endothelial
            fibroblasts
            smc
        kidney
            lymphoid
            myeloid
            tubular
        liver
            hepatocytes
            kupffer
            lymphoid
        lung
            alveolar
            macrophages
        lymph_node
            lymphoid
        ovary
            endothelial
            fibroblasts
            granulosa
            macrophages
            smc
            theca
        pancreas
            ductal
            endocrine
            exocrine-glandular
        pbmc
            lymphoid
            myeloid
        placenta
            fibroblasts
            hofbauer-cells
            trophoblasts
        prostate
            basal-prostatic
            endothelial
            fibroblasts
            prostatic-glandular
            smc
            urothelial
        rectum
            enterocytes
            goblet-cells
            paneth
        skeletal_muscle
            endothelial
            fibroblasts
            lymphoid
            myeloid
            myocytes
            smc
        skin
            endothelial
            fibroblasts
            keratinocytes
            langerhans-cells
            lymphoid
            smc
        small_intestine
            enterocytes
        spleen
            lymphoid
        stomach
            lymphoid
        testis
            endothelial
            leydig-cells
            monocyte
            peritubular-cells
            spermatids
            spermatocytes
            spermatogonia

    References
    ----------
    [1] Obradovic, A., Vlahos, L., Laise, P., Worley, J., Tan, X., Wang, A., & Califano, A. (2021). PISCES: A pipeline for the systematic, protein activity-based analysis of single cell RNA sequencing data. Biorxiv, 6, 22.
    """
    return load_pisces_network(tissue, celltype, desired_format)
