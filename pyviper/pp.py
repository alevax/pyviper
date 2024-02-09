### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
import numpy as np
import anndata
import string
from scipy.stats import rankdata
import string
import random
import scanpy as sc
from tqdm.auto import tqdm
import os

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


def __sigT(x, slope = 20, inflection = 0.5):
    return (1 - 1/(1 + np.exp(slope * (x - inflection))))

def viper_similarity(adata,
                     nn = None,
                     ws = [4, 2],
                     alternative=['two-sided','greater','less'],
                     layer=None,
                     filter_by_feature_groups=None):
    """\
    If ws is a single number, weighting is performed using an exponential function.
    If ws is a 2 numbers vector, weighting is performed with a symmetric sigmoid
    function using the first element as inflection point and the second as trend.

    Parameters
    ----------
    adata
        An anndata.AnnData containing protein activity (NES), where rows are
        observations/samples (e.g. cells or groups) and columns are features
        (e.g. proteins or pathways).
    nn (default: None)
        Optional number of top regulators to consider for computing the similarity
    ws (default: [4, 2])
        Number indicating the weighting exponent for the signature, or vector of
        2 numbers indicating the inflection point and the value corresponding to
        a weighting score of .1 for a sigmoid transformation, only used if nn is
        ommited.
    alternative (default: 'two-sided')
        Character string indicating whether the most active (greater), less
        active (less) or both tails (two.sided) of the signature should be used
        for computing the similarity.
    layer (default: None)
        The layer to use as input data to compute the signatures.
    filter_by_feature_groups (default: None)
        The selected regulators, such that all other regulators are filtered out
        from the input data. If None, all regulators will be included. Regulator
        sets must be from one of the following: "tfs", "cotfs", "sig", "surf".

    Returns
    -------
    The original anndata.AnnData object where adata.obsp['viper_similarity']
    contains a signature-based distance numpy.ndarray.

    References
    ----------
    Alvarez, M. J., Shen, Y., Giorgi, F. M., Lachmann, A., Ding, B. B., Ye, B. H.,
    & Califano, A. (2016). Functional characterization of somatic mutations in
    cancer using network-based inference of protein activity. Nature genetics,
    48(8), 838-847.
    """

    mat = _get_anndata_filtered_by_feature_group(adata, layer, filter_by_feature_groups).to_df()

    if np.min(mat)>=0 :
        mat = rankdata(mat,axis=1)
        mat = norm.ppf(mat/(np.sum(mat.isna()==False,axis = 1)+1))

    mat[mat.isna()] =0 # will this work?

    xw = mat

    if nn == None:
        if alternative == 'greater':
            xw[xw < 0] = 0
        elif alternative == 'less' :
            xw[xw > 0] = 0

        if len(ws) == 1:
            xw = np.transpose(xw)/np.max(np.abs(mat), axis = 1)
            xw = np.sign(xw) * np.abs(xw) ** ws

        else:
            ws[1] = 1/(ws[1] - ws[0]) * np.log(1/0.9 -1)
            xw = np.sign(xw) *__sigT(np.abs(mat),ws[1],ws[0]) #why it's 1, 0 instead of 0,1

    else:
        if alternative == 'greater':
            xw = rankdata(-mat,axis=1)
            mat[xw > nn] = None
        elif alternative == 'less' :
            xw = rankdata(mat,axis=1)
            mat[xw > nn] = None
        else:
            xw = rankdata(mat,axis=1)
            mat[xw > nn/2 & xw <(len(xw) - nn/2 +1)] = None

    nes = np.sqrt(np.sum(xw**2, axis = 1))
    xw = xw.transpose()/np.sum(np.abs(xw),axis = 1)

    t2 = norm.ppf(rankdata(xw.transpose(), axis = 1)/(mat.shape[1]+1))
    vp = np.matmul(t2, xw)

    vp = vp * nes

    tmp = np.array([vp.values[np.triu_indices(vp.shape[0], 1)],vp.T.values[np.triu_indices(vp.shape[0], 1)]])
    tmp = np.sum(tmp * tmp ** 2, axis=0) / np.sum(tmp ** 2, axis=0)

    vp.values[np.triu_indices(vp.shape[0], 1)] = tmp
    vp = vp.T
    vp.values[np.triu_indices(vp.shape[0], 1)] = tmp
    vp.columns = vp.index

    adata.obsp['viper_similarity'] = vp

    return adata

def aracne3_to_regulon(net_file, net_df=None, anno=None, MI_thres=0, regul_size=50, with_count_values=False, normalize_MI_per_regulon=True):
    pd.options.mode.chained_assignment = None
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

    net.sort_values(by=['regulator.values','count.values','mi.values'],ascending=[True,False,False], inplace=True)
    net = net.groupby('regulator.values').head(regul_size)

    if normalize_MI_per_regulon:
        reg_max = net.groupby(['regulator.values'])['mi.values'].transform('max')
    else:
        reg_max = net['mi.values'].max()

    # print(reg_max)
    op = pd.DataFrame({
        'regulator': net['regulator.values'],
        'target' : net['target.values'],
        'mor': net['scc.values'],
        'likelihood': net['mi.values']/reg_max
        }
    )

    return op

def select_cells(n_cells_per_metacell,probability_weights, use_decay, decay_factor, num_cells_gq, adata_gq_cells):
    ''' 
    Select a set of cells for a metacell based on specified probability weights. A good quality cell is defined as a cell with 
    silouhette score greater than 0.1.

    Parameters:
        n_cells_per_metacell (int): Number of cells to select for the metacell.
        probability_weights (numpy.ndarray): Array of probability weights for each cell.
        use_decay (bool): If True, use penalized probability sampling with decay.
        decay_factor (float): Decay factor for penalized probability.
        num_cells_gq (int): Total number of good-quality cells.
        adata_gq_cells (anndata.AnnData): Anndata object containing good-quality cells data.

    Returns:
        numpy.ndarray: Array of selected cell indices for the metacell.

    If `use_decay` is True, the function performs penalized probability sampling
    using the specified probability weights and decay factor. It selects `n_cells_per_metacell`
    cells based on a combined score of probabilities and distances. If `use_decay` is False,
    the function performs uniform probability sampling.

    The selected cell indices are returned as a numpy array.

    
    '''
    if use_decay:
        #Using Penalized Probability upon Resampling 

        #select a starting cell
        a_selected_cell = np.random.choice(num_cells_gq,p= probability_weights)
        #obtain its neighbors and the distances
        dist_array = np.ravel(adata_gq_cells.obsp['distances'].todense()[a_selected_cell, :])
        indices = np.nonzero(dist_array)[0]
        nonzero_distances = dist_array[indices]

        #extract the probabilities of those neighbors
        probabilities_cells = probability_weights[indices]


        #produce a combined score using the probabilities and the distances. 
        #Higher probability and smaller distances are prioritized.
        combined_score = probabilities_cells / nonzero_distances
        combined_score = [score / sum(combined_score) for score in combined_score]

        #select cells to create the metacell
        selected_cells_int = np.random.choice(indices, n_cells_per_metacell,p= combined_score, replace=False)
        selected_cells_int = np.append(selected_cells_int, a_selected_cell)
    else: 
        #Using Uniform Probability

        #select a starting cell
        a_selected_cell = random.sample(range(0, adata_gq_cells.shape[0]), 1)

        #obtain its neighbors and the distances
        dist_array = np.ravel(adata_gq_cells.obsp['distances'].todense()[a_selected_cell, :])

        #select the cells with the smallest distances to create the metacell
        selected_cells_int = sparse_argsort(dist_array)[0:n_cells_per_metacell] # Ascend Ordering
        selected_cells_int = np.append(selected_cells_int, a_selected_cell)
    return selected_cells_int



def sparse_argsort(array):
    indices = np.nonzero(array)[0]
    return indices[np.argsort(array[indices])]

def generateRandomCellID(prefix='',chars=string.ascii_uppercase + string.digits, N=10):
    return "CellID_" + prefix + '_' + ''.join(random.choice(chars) for _ in range(N))

def generateMetacellAnnData(arg0, n_metacells_per_cluster=500, n_metacells=1000, n_cells_per_metacell=10,
                            n_neighbors = 15,
                            cluster_label="clusters" , method="proportional", experiment_dir_path="",
                            use_decay = False, decay_factor = 0.1, seed = 42):
    '''
    Generate an AnnData object representing metacells from an input AnnData object.

    Parameters:
        arg0 (anndata.AnnData): Input AnnData object containing single-cell data.
        n_metacells_per_cluster (int): Number of metacells to generate for each cluster.
        n_metacells (int): Total number of metacells to generate.
        n_cells_per_metacell (int): Number of cells to include in each metacell.
        n_neighbors (int): Number of neighbors for computing cell-cell distances.
        cluster_label (str): Name of the column in the obs attribute specifying cell clusters.
        method (str): Method for generating metacells ("proportional", "absolute", "relative").
        experiment_dir_path (str): Path to the directory to save the generated data.
        use_decay (bool): If True, use probability decay during metacell generation.
        decay_factor (float): Decay factor for probability during metacell generation.
        seed (int): Seed for random number generation. 
    Returns:
        anndata.AnnData: AnnData object representing the generated metacells.

    Method Differences:
        - "proportional": This method generates metacells based on the proportion of cells in each cluster.
          Metacells are created proportionally to the cluster sizes, and cell probabilities may decay.
        - "absolute": Metacells are generated with an absolute number of cells per metacell.
          The method creates a fixed number of cells in each metacell, and cell probabilities may decay.
        - "relative": Metacells are generated with a relative number of cells per metacell.
          The method creates a variable number of cells based on cluster sizes, and cell probabilities may decay.

    This function generates metacells from the input single-cell data in the form of an AnnData object.
    The method of metacell generation can be selected based on the `method` parameter. The resulting
    AnnData object contains metacell data, and the generated metacells are saved in the specified
    `experiment_dir_path`.

    '''
    
    _adata = arg0.copy()
    n_var = _adata.shape[1]

    #check for missing columns
    if "silhouette_score" not in _adata.obs.columns:
        raise Exception("Please provide silhouette score as a column in adata observations.")
    gq_cell_counts = {}
    # sc.pp.neighbors(_adata, n_neighbors=51, n_pcs=30)

    var_names = _adata.var_names.tolist()

    #for the proportional approach, it is necessary to have the distribution of the clusters
    _adata.obs[cluster_label] = _adata.obs[cluster_label].astype('category')
    cluster_perc_dict = {i: [] for i in _adata.obs[cluster_label].cat.categories.to_list()}
    tot_cells = len(_adata.obs)
    for a_cluster in _adata.obs[cluster_label].cat.categories:
        index = _adata.obs[cluster_label].isin([a_cluster])
        tot_cells_per_cluster = index.sum()
        perc = round(tot_cells_per_cluster / tot_cells * 100, 1)
        cluster_perc_dict[a_cluster] = round(perc / 100 * n_metacells)
        
  

    random.seed(seed)
    list_of_matrices = list()

    for a_cluster, cluster_iterations in cluster_perc_dict.items():

        
        index = _adata.obs[cluster_label].isin([a_cluster]) & _adata.obs["silhouette_score"] > 0.1
        adata_gq_cells = _adata[index, :]

        sc.pp.neighbors(adata_gq_cells, n_neighbors=n_neighbors, n_pcs=30)
        

        num_cells_gq = adata_gq_cells.shape[0]
        for cell_id in adata_gq_cells.obs_names:
            gq_cell_counts[cell_id] = [0, a_cluster]
        print("Good quality cells ", num_cells_gq)
        #initial probabilities are discrete uniform
        probability_weights = np.full(num_cells_gq, 1/num_cells_gq)

        if method == "proportional":
            m = np.zeros([cluster_iterations, n_var])
            m = np.matrix(m)
            for a_metacell in tqdm(range(0, cluster_iterations), desc="Metacells for Cluster: " + str(a_cluster)):
                
                selected_cells_int = select_cells(n_cells_per_metacell,probability_weights,use_decay,decay_factor,num_cells_gq,adata_gq_cells)
                #once a cell is selected, decay its probability using a decay factor
                probability_weights[selected_cells_int]*=decay_factor
                #then, re-normalize the probabilities to ensure that they sum up to 1
                probability_weights = probability_weights / np.sum(probability_weights)
                selected_cell_ids = [adata_gq_cells.obs_names[i] for i in selected_cells_int]
                for selected_cell_id in selected_cell_ids:
                    gq_cell_counts[selected_cell_id][0] += 1

                tmp = adata_gq_cells.raw[adata_gq_cells.obs.index[selected_cells_int], :].X.sum(axis=0)
                m[a_metacell,] = tmp
                if (tmp.sum() == 0):
                    print("No counts found in metacells")
            list_of_matrices.append(m)

        elif method == "absolute":
            m = np.zeros([n_metacells_per_cluster, n_var])
            m = np.matrix(m)
            for a_metacell in tqdm(range(0, n_metacells_per_cluster), desc="Metacells for Cluster: " + str(a_cluster)):
                selected_cells_int = select_cells(n_cells_per_metacell,probability_weights,use_decay,decay_factor,num_cells_gq,adata_gq_cells)
                #once a cell is selected, decay its probability using a decay factor
                probability_weights[selected_cells_int]*=decay_factor
                #then, re-normalize the probabilities to ensure that they sum up to 1
                probability_weights = probability_weights / np.sum(probability_weights)
                selected_cell_ids = [adata_gq_cells.obs_names[i] for i in selected_cells_int]
                for selected_cell_id in selected_cell_ids:
                    gq_cell_counts[selected_cell_id][0] += 1
               
                tmp = adata_gq_cells.raw[adata_gq_cells.obs.index[selected_cells_int], :].X.sum(axis=0)
                m[a_metacell,] = tmp
                if (tmp.sum() == 0):
                    print("No counts found in metacells")
            list_of_matrices.append(m)
        elif method == "relative":
            #only create n_metacells_per_cluster if there are sufficient good quality cells 
            #to produce metacells without re-using too many cells
            max_limit = int(num_cells_gq/n_cells_per_metacell)
            max_limit = min(max_limit,n_metacells_per_cluster)
            m = np.zeros([max_limit, n_var])
            m = np.matrix(m)
            for a_metacell in tqdm(range(0, max_limit), desc="Metacells for Cluster: " + str(a_cluster)):
                selected_cells_int = select_cells(n_cells_per_metacell,probability_weights,use_decay,decay_factor,num_cells_gq,adata_gq_cells)
                #once a cell is selected, decay its probability using a decay factor
                probability_weights[selected_cells_int]*=decay_factor
                #then, re-normalize the probabilities to ensure that they sum up to 1
                probability_weights = probability_weights / np.sum(probability_weights)
                selected_cell_ids = [adata_gq_cells.obs_names[i] for i in selected_cells_int]
                for selected_cell_id in selected_cell_ids:
                    gq_cell_counts[selected_cell_id][0] += 1
                
                tmp = adata_gq_cells.raw[adata_gq_cells.obs.index[selected_cells_int], :].X.sum(axis=0)
                m[a_metacell,] = tmp
                if (tmp.sum() == 0):
                    print("No counts found in metacells")
            list_of_matrices.append(m)
        else:
            print("ERROR")
        [print(mat.shape) for mat in list_of_matrices]


    i = 0
    cluster_labels = _adata.obs[cluster_label].cat.categories.to_list()
    obs_names_full = []
    folder = "with_penalty" if use_decay else "without_penalty"
    output_path = os.path.join(experiment_dir_path, folder, str(n_cells_per_metacell))
    if not os.path.exists(output_path):
        os.makedirs(output_path) 

    for mat in list_of_matrices:
        sample_name = "C_" + str(i)
        cluster_name = "C_" + str(cluster_labels[i])
        cluster_name = cluster_name.replace(' ', '_')
        
        i = i+1
        obs_names = [generateRandomCellID(prefix=sample_name,N=10) for i in range(0, mat.shape[0])]
        obs_names = obs_names[0:mat.shape[0]]
        obs_names_full.append(obs_names)
        x = anndata.AnnData(X=mat.A, obs=pd.DataFrame(index=obs_names), var=pd.DataFrame(index=var_names))
        sc.pp.normalize_total(x, inplace=True, target_sum=1e6)
        #x.write(filename)
        df = pd.DataFrame(data=x.X.T, index=var_names, columns=obs_names)
        filename = os.path.join(output_path ,"whole-dataset-metacells-expression_"+ cluster_name + "-expression.txt")
        df.insert(0, "Gene", var_names)
        df.to_csv(filename, sep = "\t", index = False)

    m = np.vstack(list_of_matrices)
    repeated_cluster_labels = [label for label, mat in zip(cluster_labels, list_of_matrices) for _ in range(mat.shape[0])]

    print(m.shape)


    obs_names_full_flat = [item for sublist in obs_names_full for item in sublist]
    adata_metacell = anndata.AnnData(X=m.A, obs=pd.DataFrame(index=obs_names_full_flat), var=pd.DataFrame(index=var_names))
    adata_metacell.obs[cluster_label] = repeated_cluster_labels
    adata_metacell.raw = adata_metacell
    sc.pp.normalize_total(adata_metacell, inplace=True, target_sum=1e6)

    # filename = os.path.join(experiment_dir_path, sample_name + "-whole-dataset-metacells-expression.h5ad")

    filename = os.path.join(output_path ,"whole-dataset-metacells-expression_" + str(n_cells_per_metacell) + ".h5ad")
    pd.DataFrame.from_dict(gq_cell_counts, orient='index', columns=['Count', 'Cluster']).to_csv(os.path.join(output_path,"metacell_counts" + "_" + str(n_cells_per_metacell) + ".csv"))
    adata_metacell.write(filename)
    return (adata_metacell)
