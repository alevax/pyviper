### ---------- IMPORT DEPENDENCIES ----------
import pandas as pd
from scipy.stats import rankdata
from scipy.stats import norm
import numpy as np
import warnings

### ---------- EXPORT LIST ----------
__all__ = ['NaRnEA_classic']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------ HELPER FUNCTIONS -----------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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
    # a = np.sqrt(- final_p + np.log(2))
    # pos_nes = -1.4374174 + 1.8396835*a - 0.0562393*a**2 + 0.0025810*a**3
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

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# -----------------------------------------------------------------------------
# ------------------------------- MAIN FUNCTION -------------------------------
# -----------------------------------------------------------------------------
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def NaRnEA_classic(gex_data,
                   interactome,
                   layer=None,
                   eset_filter=False,
                   min_targets=30,
                   verbose=False,
                   return_as_df=False):

    # filter out those with target less than min.targets
    interactome = interactome.copy()
    interactome.prune(min_targets = min_targets, max_targets = None, eliminate = False, verbose = False)

    if (eset_filter):
        # This will affect the rankings of genes by eliminating those not present in the interactome
        tmp = np.unique(np.concatenate((interactome.get_targetSet(), interactome.get_regulonNames())))
        gex_data = gex_data[:,gex_data.var_names.isin(pd.Series(tmp))]

    pd.options.mode.chained_assignment = None
    exp_genes = list(gex_data.var.sort_index().index)


    # modify regulon list, take the intersect of targets and genes, making there is no 0 nor +_1 in am
    # filtered_table = int_table[int_table['target'].isin(exp_genes)]
    n_targets_not_in_exp_genes = np.count_nonzero(~np.isin(interactome.get_targetSet(), exp_genes))
    if n_targets_not_in_exp_genes > 0:
        # raise ValueError(
        warnings.warn('interactome "' + str(interactome.name) + '" contains ' +
                         str(n_targets_not_in_exp_genes) + " targets missing from gex_data.var.\n\t" +
                        "Please run interactome.filter_targets(gex_data.var_names) on your network to\n\t" +
                         "resolve this. It is highly recommend to do this on the unPruned network and\n\t"+
                         "then prune to the pruned network contains a consistent number of targets per\n\t"
                         "regulator, allow of which exist within gex_data.")
        interactome.filter_targets(gex_data.var_names)
    int_table = interactome.net_table
    filtered_table = int_table

    # why I need to remove that, there is no filtered_table[filtered_table['mor'] == -1]
    filtered_table['mor'].replace(1 ,0.999, inplace= True)
    filtered_table['mor'].replace(-1 ,-0.999, inplace= True)
    filtered_table['mor'] = filtered_table['mor'].apply(lambda x: replace_random(x, -0.001, 0.001))

    #filtered_table['mor'].replace(0 ,0.999, inplace= True)
    # why we cant have 0, 1 and -1?

    # modify the expression matrix:
        # fill with number less than the min in this sample, sign is randomly assigned according to the +_ proportion
    if verbose: print('reordering genes')

    if layer is None:
        expr = np.copy(gex_data[:,exp_genes].X)
    else:
        expr = np.copy(gex_data[:,exp_genes].layers[layer])



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
    genes = list(gex_data.var.index)
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

    if verbose: print('Calculating DES...')
    D_list = directed_nes(ranks, signs, AW_AM_prob, E_rs, E_r2, expr.shape[0])

    if verbose: print('Calculating UES...')
    U_list = undirected_nes(ranks, signs, AW_AM_abs_prob, E_r, E_r2, expr.shape[0])

    if verbose: print('Calculating NES...')
    COV_nes = nes_covariance(ranks, signs, AW_AM_prob, AW_AM_abs_prob, E_r, E_rs,D_list[2], U_list[2], len(genes))
    NES_mat = combine_nes(D_list[3], U_list[3], COV_nes)

    if verbose: print('Calculating PES...')
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



    if (return_as_df) == False :
        NES_mat = pd.DataFrame(NES_mat.T, index = gex_data.obs.index, columns = list(AW_AM_prob.columns))
        PES_mat = pd.DataFrame(PES_mat.T, index = gex_data.obs.index, columns = list(AW_AM_prob.columns))


    result = {"nes": NES_mat, "pes": PES_mat}
    return result
