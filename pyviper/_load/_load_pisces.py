### Import dependencies
import pandas as pd
from ..interactome import Interactome

### EXPORT LIST
__all__ = ['load_pisces_network']

pisces_networks_dict = {
    'adipose_tissue' : [
        'adipocytes',
        'dendritic',
        'fibroblasts',
        'lymphoid',
        'myeloid',
        'smc'
    ],
    'bone_marrow' : ['lymphoid'],
    'brain' : ['neuron'],
    'breast': [
        'adipocytes',
        'dendritic',
        'endothelial',
        'epithelial',
        'fibroblasts',
        'lymphoid',
        'myeloid',
        'smc'
    ],
    'bronchus' : ['epithelial', 'keratinocytes'],
    'colon' : ['enterocyte', 'goblet', 'lymphoid'],
    'endometrium' : [
        'adipocytes',
        'endothelial',
        'epithelial',
        'fibroblasts',
        'lymphoid'
    ],
    'esophagus': ['epithelial', 'fibroblasts'],
    'eye' : ['bipolar', 'muller-glia', 'photoreceptor'],
    'heart_muscle' : ['cardiomyocyte', 'endothelial', 'fibroblasts', 'smc'],
    'kidney' : ['lymphoid', 'myeloid', 'tubular'],
    'liver' : ['hepatocytes', 'kupffer', 'lymphoid'],
    'lung' : ['alveolar', 'macrophages'],
    'lymph_node': ['lymphoid'],
    'ovary' : ['endothelial', 'fibroblasts', 'granulosa', 'macrophages', 'smc', 'theca'],
    'pancreas' : ['ductal', 'endocrine', 'exocrine-glandular'],
    'pbmc' : ['lymphoid', 'myeloid'],
    'placenta' : [
        'fibroblasts', 'hofbauer-cells', 'trophoblasts'
    ],
    'prostate' : [
        'basal-prostatic',
        'endothelial',
        'fibroblasts',
        'prostatic-glandular'
        'smc'
        'urothelial'
    ],
    'rectum' : ['enterocytes', 'goblet-cells', 'paneth'],
    'skeletal_muscle' : [
        'endothelial',
        'fibroblasts',
        'lymphoid',
        'myeloid',
        'myocytes',
        'smc'
    ],
    'skin' : [
        'endothelial',
        'fibroblasts',
        'keratinocytes',
        'langerhans-cells',
        'lymphoid',
        'smc'
    ],
    'small_intestine' : ['enterocytes'],
    'spleen' : ['lymphoid'],
    'stomach' : ['lymphoid'],
    'testis' : [
        'endothelial',
        'leydig-cells',
        'monocyte',
        'peritubular-cells',
        'spermatids',
        'spermatocytes',
        'spermatogonia'
    ]
}

def load_parquet_from_github(url):
    """
    Load a parquet file from a public GitHub repository.

    Args:
        url (str): The URL of the parquet or parquet.gzip file on GitHub.

    Returns:
        pandas.DataFrame: The loaded parquet data as a DataFrame.
    """
    try:
        # Read the CSV file from the provided URL
        df = pd.read_parquet(url)
        return df
    except Exception as e:
        print("An error occurred:", e)
        return None

def load_pisces_network(tissue, celltype, desired_format = "human_ensembl"):
    available_tissues = list(pisces_networks_dict.keys())
    if tissue not in available_tissues:
        raise ValueError('tissue "' + str(tissue) + '" not available. Choose from:\n\t' + "\n\t".join(available_tissues) )

    available_celltypes = pisces_networks_dict.get(tissue)

    if celltype not in available_celltypes:
        raise ValueError('celltype "' + str(celltype) + '" of tissue "' + str(tissue) + '" not available. Choose from:\n\t' + "\n\t".join(available_celltypes) )

    url = 'https://raw.githubusercontent.com/califano-lab/PISCES_networks/main/parquet/'
    url += tissue + "_" + celltype + '.parquet.gzip'

    net_table = load_parquet_from_github(url)

    net = Interactome(tissue + "_" + celltype, net_table)

    if desired_format != "human_ensembl":
        net.translate_regulators(desired_format, verbose = False)
        net.translate_targets(desired_format, verbose = False)

    return net
