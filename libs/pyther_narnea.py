### Import  dependencies

import pandas as pd
from pyther_classes import *
from scipy.stats import norm
from scipy.stats import rankdata
from joblib import Parallel, delayed
import anndata
### Tool functions


# def log_quantile(log_r):
#     A0 = 3.3871327179e0
#     A1 = 5.0434271938e1
#     A2 = 1.5929113202e2
#     A3 = 5.9109374720e1
#     B1 = 1.7895169469e1
#     B2 = 7.8757757664e1
#     B3 = 6.7187563600e1
#
#     C0 = 1.4234372777e0
#     C1 = 2.7568153900e0
#     C2 = 1.3067284816e0
#     C3 = 1.7023821103e-1
#     D1 = 7.3700164250e-1
#     D2 = 1.2021132975e-1
#
#     E0 = 6.6579051150e0
#     E1 = 3.0812263860e0
#     E2 = 4.2868294337e-1
#     E3 = 1.7337203997e-2
#     F1 = 2.4197894225e-1
#     F2 = 1.2258202635e-2
#
#     SPLIT1 = 0.425
#     CONST1 = 0.180625
#     SPLIT2 = 5.0e0
#     CONST2 = 1.6e0
#
#     if log_r > np.log(0.075):
#         p = np.exp(log_r)
#         q = p - 0.5
#         R = CONST1 - q * q
#         PPND7 = -q * (((A3 * R + A2) * R + A1) * R + A0) / (((B3 * R + B2) * R + B1) * R + 1)
#     else:
#         R = np.sqrt(-log_r)
#         if R <= SPLIT2:
#             R = R - CONST2
#             PPND7 = (((C3 * R + C2) * R + C1) * R + C0) / ((D2 * R + D1) * R + 1)
#         else:
#             R = R - SPLIT2
#             PPND7 = (((E3 * R + E2) * R + E1) * R + E0) / ((F2 * R + F1) * R + 1)
#
#     return PPND7

def vect_quantile_fast(log_r):
    A0 = 3.3871327179e0
    A1 = 5.0434271938e1
    A2 = 1.5929113202e2
    A3 = 5.9109374720e1
    B1 = 1.7895169469e1
    B2 = 7.8757757664e1
    B3 = 6.7187563600e1

    C0 = 1.4234372777e0
    C1 = 2.7568153900e0
    C2 = 1.3067284816e0
    C3 = 1.7023821103e-1
    D1 = 7.3700164250e-1
    D2 = 1.2021132975e-1

    E0 = 6.6579051150e0
    E1 = 3.0812263860e0
    E2 = 4.2868294337e-1
    E3 = 1.7337203997e-2
    F1 = 2.4197894225e-1
    F2 = 1.2258202635e-2

    # SPLIT1 = 0.425
    CONST1 = 0.180625
    SPLIT2 = 5.0e0
    CONST2 = 1.6e0

    log_r = log_r.astype(float)

    cond_1 = log_r > np.log(0.075)
    # cond_2 = (np.sqrt(-log_r) <= SPLIT2)

    # Initialize cond_2 as an array of True values with the same shape as log_r
    cond_2 = np.ones_like(log_r, dtype=bool)
    # Calculate cond_2 only where cond_1 is not true to prevent any warnings
    cond_2[~cond_1] = (np.sqrt(-log_r[~cond_1]) <= SPLIT2)

    coords_1 = np.where(cond_1)
    coords_2 = np.where(~cond_1 & cond_2)
    coords_3 = np.where(~cond_1 & ~cond_2)

    # if log_r > np.log(0.075):
    p = np.exp(log_r[coords_1])
    q = p - 0.5
    R = CONST1 - q * q
    PPND7 = -q * (((A3 * R + A2) * R + A1) * R + A0) / (((B3 * R + B2) * R + B1) * R + 1)
    log_r[coords_1] = PPND7
    # elif R <= SPLIT2:
    R = np.sqrt(-log_r[coords_2])
    R = R - CONST2
    PPND7 = (((C3 * R + C2) * R + C1) * R + C0) / ((D2 * R + D1) * R + 1)
    log_r[coords_2] = PPND7
    # else:
    R = np.sqrt(-log_r[coords_3])
    R = R - SPLIT2
    PPND7 = (((E3 * R + E2) * R + E1) * R + E0) / ((F2 * R + F1) * R + 1)
    log_r[coords_3] = PPND7

    return log_r

def directed_nes(rank_mat, sign_mat, AW_AM_prob, E_rs, E_r2, nsample):
    D = np.matmul(rank_mat*sign_mat, AW_AM_prob) #samples X regulators
    D_e = np.matmul(np.matrix(np.sum(AW_AM_prob, axis = 0)).transpose(),np.matrix(E_rs)) #621(regs)*211(samples)
    reg_squared_sum_product = np.matrix(np.sum(AW_AM_prob**2, axis = 0)) #1*612
    E_d_2 = reg_squared_sum_product.transpose() @ np.matrix(E_rs**2)
    E_d2 =  reg_squared_sum_product.transpose() @ np.matrix(np.repeat(E_r2, nsample))
    D_V = E_d2 - E_d_2
    #print(D.shape)
    #print(D_e.shape)
    #print(D_V.shape)

    D_nes = (D.transpose() - D_e)/np.sqrt(D_V)
    return([D, D_e, D_V, D_nes])


def undirected_nes(rank_mat, sign_mat, AW_AM_abs_prob, E_r, E_r2, nsample):
    U = np.matmul(rank_mat, AW_AM_abs_prob)
    U_e = np.matmul(np.matrix(np.sum(AW_AM_abs_prob, axis = 0)).transpose(),np.matrix(E_r)) #621(regs)*211(samples)
    reg_squared_sum_product_abs = np.matrix(np.sum(AW_AM_abs_prob**2, axis = 0)) #1*612
    E_u_2 =  reg_squared_sum_product_abs.transpose() @ np.matrix(E_r**2)
    E_u2 =  reg_squared_sum_product_abs.transpose() @ np.matrix(np.repeat(E_r2, nsample))
    U_V = E_u2 - E_u_2
    U_nes = (U.transpose() - U_e)/np.sqrt(U_V)
    return([U,U_e, U_V, U_nes])


def nes_covariance(rank_mat, sign_mat,AW_AM_prob, AW_AM_abs_prob, E_r, E_rs, D_V, U_V, ngene):
    cov_prod_mat = np.matrix(np.sum(AW_AM_prob * AW_AM_abs_prob,axis = 0))
    E_r2s =  np.matrix((1/ngene) * np.sum( rank_mat * rank_mat * sign_mat, axis= 1))
    E_du = cov_prod_mat.transpose() @ E_r2s
    E_d_u = cov_prod_mat.transpose() @ np.matrix(E_rs * E_r)
    COV_du = E_du - E_d_u

    COV_nes = COV_du/np.sqrt(np.multiply(D_V , U_V))
    return(COV_nes)


def combine_nes(D_nes, U_nes, COV_nes):

    NES_pos = (D_nes + U_nes) / np.sqrt(2 + 2 * COV_nes)
    NES_neg = (D_nes - U_nes) / np.sqrt(2 - 2 * COV_nes)


    p_pos = norm.logsf(NES_pos)
    p_neg = norm.logcdf(NES_neg)
    p_dif = (p_pos < p_neg)
    min_p = (p_pos * p_dif + p_neg * (~p_dif))
    final_p = min_p + np.log(2) + np.log1p(np.exp(min_p) / (-2))

    # facing precision issue with this approach, switch to a numerical approximation way

    # quantile =  np.exp(final_p - np.log(2))
    # pos_nes = norm.isf(quantile, loc=0, scale=1)
    # neg_nes = norm.ppf(quantile, loc=0, scale=1)

    # a = np.sqrt(- final_p + np.log(2))
    # pos_nes = -1.4374174 + 1.8396835*a - 0.0562393*a**2 + 0.0025810*a**3

    # a = final_p - np.log(2)
    # vect_quantile = np.vectorize(log_quantile)
    # vect_quantile = np.vectorize(log_quantile_fast)
    # pos_nes = vect_quantile(a)

    a = final_p - np.log(2)
    pos_nes = vect_quantile_fast(a)
    neg_nes = -pos_nes

    NES_mat = (pos_nes * p_dif + neg_nes * (~p_dif))
    return NES_mat


def replace_random(x, a, b):
    if x == 0:
        rand = np.random.random()# 0,1

        return rand*(b-a) + a
    else:
        return x


### Matrix narnea


def matrix_narnea(gesObj, intObj, intermediate = False, min_targets = 30,verbose = False):

    int_table = intObj.netTable

    pd.options.mode.chained_assignment = None
    exp_genes = list(gesObj.var.sort_index().index)


    # modify regulon list, take the intersect of targets and genes, making there is no 0 nor +_1 in am

    filtered_table = int_table[int_table['target'].isin(exp_genes)]
    # why I need to remove that, there is no filtered_table[filtered_table['mor'] == -1]
    filtered_table['mor'].replace(1 ,0.999, inplace= True)
    filtered_table['mor'].replace(-1 ,-0.999, inplace= True)
    filtered_table['mor'] = filtered_table['mor'].apply(lambda x: replace_random(x, -0.001, 0.001))
    #filtered_table['mor'].replace(0 ,0.999, inplace= True)
    # why we cant have 0, 1 and -1?

    # filtered out those with target less than min.targets )may be this should  be a function for interactone class
    reg_counts  = filtered_table['regulator'].value_counts()
    reg_counts = list(reg_counts[reg_counts>=min_targets].index)
    filtered_table = filtered_table[filtered_table['regulator'].isin(reg_counts)]

    # modify the expression matrix:
        # fill with number less than the min in this sample, sign is randomly assigned according to the +_ proportion
    if verbose:
        print('reordering genes')

    expr = np.copy(gesObj[:,exp_genes].X)


    zeros = len(expr[expr == 0])

    pos = np.sum(np.sign(expr) + np.abs(np.sign(expr)))/2
    total = np.sum(np.abs(np.sign(expr)))
    ratio = pos/total

    posnum = int(ratio * zeros )
    negnum = zeros - posnum

    sign = np.hstack((np.repeat(-1,negnum),np.repeat(1,posnum)))
    np.random.shuffle(sign)
    value = np.random.rand(zeros) * np.min(expr[expr != 0])


    expr[expr == 0] = sign * value

    # rank the gene expression in each sample and preserve the sign

    ranks = rankdata(np.abs(expr), axis = 1)
    signs = np.sign(expr)


    # initializing matrixs with
    genes = list(gesObj.var.index)
    regulons = list(filtered_table['regulator'].drop_duplicates())

    # wait, not necesscarily do that, can we do sth like merge and pivot longer?
    dat1 = pd.DataFrame({'target':genes})

    dat2 = dat1.merge(filtered_table, on='target', how = 'left')
    # should I use left or inner? #no. have to use left join
    #dat2.dropna(inplace= True)
    AM_mat = dat2.pivot(index='target',columns='regulator',values = 'mor').fillna(0)
    AM_mat = AM_mat[AM_mat.columns.dropna()]


    AW_mat =  dat2.pivot(index='target',columns='regulator',values = 'likelihood').fillna(0)
    AW_mat = AW_mat[AW_mat.columns.dropna()]

    # fill the matrix with genes == targets
    E_r = (len(genes) + 1)/2
    E_r2 = (2 * len(genes)**2 + 3 * len(genes) + 1)/6
    E_rs = 1/len(genes) * np.sum(ranks * signs, axis = 1)

    AM_abs_mat = 1-np.abs(AM_mat)
    AM_abs_mat[AM_mat == 0] = 0 #why it's like this ？

    AW_AM_prob = AM_mat * AW_mat
    AW_AM_abs_prob = AM_abs_mat * AW_mat

#    AW_AM_prob.astype('float64')
#    AW_AM_abs_prob.astype('float64')

    if verbose:
        print('Calculating DES...')

    D_list = directed_nes(ranks, signs, AW_AM_prob, E_rs, E_r2, expr.shape[0])

    if verbose:
        print('Calculating UES...')

    U_list = undirected_nes(ranks, signs, AW_AM_abs_prob, E_r, E_r2, expr.shape[0])

    COV_nes = nes_covariance(ranks, signs, AW_AM_prob, AW_AM_abs_prob, E_r, E_rs,D_list[2], U_list[2], len(genes))

    if verbose:
        print('Calculating NES...')

    NES_mat = combine_nes(D_list[3], U_list[3], COV_nes)

    if verbose:
        print('Calculating PES...')


    # for each regulator, calculate two values

    filtered_table['rank'] = filtered_table.groupby('regulator')['likelihood'].rank(method = 'dense', ascending = False)
    filtered_table['rank'] = len(genes) - filtered_table['rank'] + 1

    filtered_table['abs_mor'] = np.abs(filtered_table['mor'])

    filtered_table['d'] = filtered_table['abs_mor']* filtered_table['likelihood'] * filtered_table['rank']
    filtered_table['u'] = (1- filtered_table['abs_mor'])* filtered_table['likelihood'] * filtered_table['rank']

    #max_du = pd.concat([filtered_table.groupby(['regulator'])['d'].sum(), filtered_table.groupby(['regulator'])['u'].sum()], axis=1)
    #max_du
    # if I want it I can directly create the two matrix
    PES_pos_D  = filtered_table.groupby(['regulator'])['d'].sum().values
    PES_pos_D = PES_pos_D[:, np.newaxis] * np.ones(expr.shape[0])[np.newaxis, :]


    PES_pos_U = filtered_table.groupby(['regulator'])['u'].sum().values
    PES_pos_U = PES_pos_U[:, np.newaxis] * np.ones(expr.shape[0])[np.newaxis, :]

    PES_pos_D_nes = (PES_pos_D - D_list[1])/np.sqrt(D_list[2])
    PES_pos_U_nes = (PES_pos_U - U_list[1])/np.sqrt(U_list[2])
    PES_pos_NES = combine_nes(PES_pos_D_nes, PES_pos_U_nes, COV_nes)

    PES_neg_D = PES_pos_D * (-1)
    PES_neg_U = PES_pos_U
    PES_neg_D_nes = (PES_neg_D - D_list[1])/np.sqrt(D_list[2])
    PES_neg_U_nes = (PES_neg_U - U_list[1])/np.sqrt(U_list[2])
    PES_neg_NES = combine_nes(PES_neg_D_nes, PES_neg_U_nes, COV_nes)

    pos_NES = (NES_mat > 0)
    PES_comb_nes = PES_pos_NES * pos_NES + PES_neg_NES * (~pos_NES)
    PES_mat = NES_mat / np.abs(PES_comb_nes)



    if (intermediate) == False :
        NES_mat = pd.DataFrame(NES_mat.T, index = gesObj.obs.index, columns = list(AW_AM_prob.columns))
        PES_mat = pd.DataFrame(PES_mat.T, index = gesObj.obs.index, columns = list(AW_AM_prob.columns))



    result = {"nes": NES_mat, "pes": PES_mat}
    return result


def get_pes_list(results):
    # Iterate through each DataFrame to collect unique gene names
    pes_list = []
    for res in results:
        pes_list.append(res['pes'])
    return pes_list

def get_all_regs_in_NaRnEA_list(results):
    pes_list = get_pes_list(results)
    return get_all_regs(pes_list)

def get_all_regs(mat_list):
    # Iterate through each DataFrame to collect unique gene names
    all_regs_set = set()
    for mat in mat_list:
        reg_names = mat.columns.tolist()
        all_regs_set.update(reg_names)

    # Sort the gene names alphabetically
    all_regs = sorted(list(all_regs_set))

    return all_regs

def get_resized_mats(mat_list, empty_value = np.nan):
    # Get names of all regulators
    all_regs = get_all_regs(mat_list)

    # Iterate through each DataFrame to normalize it
    resized_mat_list = []
    for mat in mat_list:
        # Create a DataFrame with empty values (e.g. NaN or 0) for missing gene names
        mat = mat.copy()
        missing_regs = list(set(all_regs) - set(mat.columns))
        empty_df = pd.DataFrame(empty_value, index=mat.index, columns=list(missing_regs))
        # Concatenate the original DataFrame with the empty (e.g. NaN or 0) DataFrame
        resized_mat = pd.concat([mat, empty_df], axis=1)
        # Sort the columns alphabetically
        resized_mat = resized_mat[all_regs]
        resized_mat_list.append(resized_mat)

    return resized_mat_list

def get_resized_pes(results):
    # Get list of PES dataframes
    pes_list = get_pes_list(results)
    resized_pes_list = get_resized_mats(pes_list, empty_value = np.nan)
    return resized_pes_list

def find_max_absolute_value_df(dataframes):
    # Stack the DataFrames into a 3D NumPy array
    stacked_arrays = [df.values for df in dataframes]
    stacked_data = np.stack(stacked_arrays, axis=-1)

    # Calculate the absolute maximum along the last axis (axis=-1)
    absolute_max_indices = np.nanargmax(np.abs(stacked_data), axis=-1).astype(float)

    # Calculate the existance of any nan values along the last axis (axis=-1)
    nan_present_indices = np.isnan(stacked_data).any(axis=-1)

    # Replace values in absolute_max_indices with NaN wherever NaN is present in nan_present_indices
    absolute_max_indices[nan_present_indices] = np.nan

    # Create a DataFrame with the same row and column names
    result_df = pd.DataFrame(absolute_max_indices, index=dataframes[0].index, columns=dataframes[0].columns)

    return result_df

def calculate_value_proportions(result_df):
    # Convert the DataFrame to a NumPy array
    result_array = result_df.values

    # Calculate the maximum value in the result array, ignore NaN values
    max_value = int(np.nanmax(result_array))

    # Initialize an empty NumPy array to store proportions
    proportions_array = np.zeros((result_array.shape[0], max_value + 1))

    # Count the occurrences of each value in each row, ignore NaN values
    for i, row in enumerate(result_array):
        integer_values = row[np.isfinite(row) & (row == np.floor(row))]
        unique_values, counts = np.unique(integer_values, return_counts=True)
        proportions_array[i, unique_values.astype(int)] = counts / len(integer_values)

    # Convert the proportions array to a DataFrame
    proportions_df = pd.DataFrame(proportions_array, columns=range(max_value + 1), index=result_df.index)

    return proportions_df

def get_net_weight(results):
    # Resize the PES matrices so they have the same size, column names and row names
    resized_pes_list = get_resized_pes(results)
    # Stack the resized PES matrices: along the z axis, calculate position of max abs value
    max_abs_vals_df = find_max_absolute_value_df(resized_pes_list)
    # For each sample, calculate the proportion of max abs values from each network
    net_weight = calculate_value_proportions(max_abs_vals_df).sort_index()

    net_weight.index.name = 'index'
    net_weight.columns.name = 'net'

    return net_weight

def integrate_NaRnEA_xes_mats(results, bg_matrix, net_weight, xes_type = 'pes'):

    pre_xes = bg_matrix.copy()
    xes_dom = bg_matrix.copy()
    n_regs = bg_matrix.shape[1]

    for i in range(0,len(results)):
        all_regs_xes = bg_matrix.copy() + results[i][xes_type]
        all_regs_xes.fillna(0, inplace=True)
        all_regs_xes.sort_index(inplace=True)
        if xes_type == 'nes':
            xes_dom = xes_dom + ( (all_regs_xes!=0) * net_weight[i].to_numpy()[:, np.newaxis] * np.ones((1, n_regs)) )**2
        else: #xes_type == 'pes':
            xes_dom = xes_dom + ( (all_regs_xes!=0) * net_weight[i].to_numpy()[:, np.newaxis] * np.ones((1, n_regs)) )
        net_weight_array = net_weight[i].to_numpy()[:, np.newaxis] * np.ones((1, n_regs))
        pre_xes =  pre_xes + all_regs_xes * net_weight_array

    if xes_type == 'nes':
        pre_xes = pre_xes/np.sqrt(xes_dom)
    else: #xes_type == 'pes':
        pre_xes = pre_xes/xes_dom

    # There must be regulator activity in some cells where neither network provides nonzero value.
    # Imagine we have a regulator A. A has a regulon A_1 and A_2 in net_1 and net_2 respectively.
    # Imagine, the targets of A_1 and the targets of A_2 are both 0 in samples S1, S2, ...., Sn.
    # For these samples pre_xes will have 0 and these 0s will be present in xes_dom in the same place.
    # This is because we start with bg_matrix (0 matrix) and add to it wherever the xes matrices are not 0 (all_regs_xes!= 0)
    # so if the xes matrices have 0s in them, this will result in 0s remaining in the dom_matrix.
    # Hence fillna fixes 0/0, which should just remain at 0.
    # The following print statements will show the same numbers of 0s.
    # print(np.count_nonzero(pre_xes == 0))
    # print(np.count_nonzero(xes_dom == 0))
    xes = pre_xes.fillna(0)

    return(xes)


### Meta narnea
def meta_narnea(gesObj, intObj, sample_weight = True, njobs = 1, verbose = True):
    pd.options.mode.chained_assignment = None

    if type(intObj) == Interactome:
        return matrix_narnea(gesObj, intObj)
    # elif len(intObj) == 1:
    #     return matrix_narnea(gesObj, intObj[0])
    elif njobs == 1:
        # results = [matrix_narnea(gesObj, iObj, verbose = verbose) for iObj in intObj]
        results = []
        tot_nets = len(intObj)
        n_completed_nets = 0
        if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")
        for iObj in intObj:
          results.append(matrix_narnea(gesObj, iObj, verbose = verbose))
          n_completed_nets = n_completed_nets + 1
          if verbose: print(str(n_completed_nets) + "/" + str(tot_nets) + " networks complete.")

    else:
        results = Parallel(n_jobs = njobs)(
            (delayed)(matrix_narnea)(gesObj,iObj)
            for iObj in intObj
            )

    if verbose:
        print('Integrating results')

    net_weight = get_net_weight(results)
    all_regs = get_all_regs_in_NaRnEA_list(results)
    bg_matrix = pd.DataFrame(0,index = results[0]['nes'].index, columns = all_regs)

    nes = integrate_NaRnEA_xes_mats(results, bg_matrix, net_weight, xes_type = 'nes')
    pes = integrate_NaRnEA_xes_mats(results, bg_matrix, net_weight, xes_type = 'pes')

    result = {"nes": nes, "pes": pes}
    return result
