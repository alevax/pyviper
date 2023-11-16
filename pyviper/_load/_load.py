### Import dependencies
from ._load_translate import *
from ._load_regulators import *
from ._load_msigdb import *

### EXPORT LIST
__all__ = [#'mouse2human',
           'human2mouse',
           'TFs', 'coTFs', 'sig', 'surf',
           'msigdb_regulon']

# def mouse2human():
#     """\
#     Retrieves the mouse to human translation pd.DataFrame from pyVIPER's data
#     folder. This dataframe contains four columns:
#         human_symbol, human_ensembl, mouse_ensembl, mouse_symbol
#
#     Returns
#     -------
#     A dataframe of :class:`~pandas.core.frame.DataFrame`.
#     """
#     return load_mouse2human()

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
    species (default: None)
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_tfs (default: None)
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
    species (default: None)
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_cotfs (default: None)
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
    species (default: None)
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_sig (default: None)
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
    species (default: None)
        When left as None, the species setting in pyviper.config will be used.
        Otherwise, manually specify "human" or "mouse".
    path_to_sig (default: None)
        When left as None, the path to surf setting in pyviper.config will be used.
        Otherwise, manually specify a filepath to a .txt file containing surface
        proteins, one on each line.

    Returns
    -------
    A list containing signaling proteins.
    """
    return load_surf(species, path_to_surf)

def msigdb_regulon(collection):
    """\
    Retrieves an object or a list of objects of class Interactome from pyviper's
    data folder containing a set of pathways from the Molecular Signatures
    Database (MSigDB), downloaded from https://www.gsea-msigdb.org/gsea/msigdb.
    These collections can be from one of the following:
        'h' for Hallmark gene sets. Coherently expressed signatures derived by
        aggregating many MSigDB gene sets to represent well-defined biological
        states or processes.
        'c2' for curated gene sets. From online pathway databases, publications
        in PubMed, and knowledge of domain experts.
        'c5' for ontology gene sets. Consists of genes annotated by the same
        ontology term.
        'c6' for oncogenic signature gene sets. Defined directly from microarray
        gene expression data from cancer gene perturbations.
        'c7' for immunologic signature gene sets. Represents cell states and
        perturbations within the immune system.

    Parameters
    ----------
    collection
        A individual string or a list of strings containing the following:
        ["h", "c2", "c5", "c6", "c7"]
        corresponding to the collections above.

    Returns
    -------
    An individual object or list of objects of class pyviper.interactome.Interactome.
    """
    return load_msigdb_regulon(collection)
