### Import  dependencies

import pandas as pd
from pyther_classes import *
from scipy.stats import norm
from scipy.stats import rankdata

### Tool functions


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
    a = np.sqrt(- final_p + np.log(2))
    pos_nes = -0.2136 + 1.417*a
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

def matrix_narnea(gesObj, int_table, intermediate = False):
    
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


    # modify the expression matrix: 
        # fill with number less than the min in this sample, sign is randomly assigned according to the +_ proportion

    expr = np.copy(gesObj[:,exp_genes].X)

    print('reordered genes')

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

    print('Calculating DES...')

    D_list = directed_nes(ranks, signs, AW_AM_prob, E_rs, E_r2, expr.shape[0])

    print('Calculating UES...')

    U_list = undirected_nes(ranks, signs, AW_AM_abs_prob, E_r, E_r2, expr.shape[0])

    COV_nes = nes_covariance(ranks, signs, AW_AM_prob, AW_AM_abs_prob, E_r, E_rs,D_list[2], U_list[2], len(genes))

    print('Calculating NES...')

    NES_mat = combine_nes(D_list[3], U_list[3], COV_nes)

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
        NES_mat = pd.DataFrame(NES_mat.T, index = gesObj.obs.index, columns= list(AW_AM_prob.columns))
        PES_mat = pd.DataFrame(PES_mat.T, index = gesObj.obs.index, columns= list(AW_AM_prob.columns))
                               

    
    result = [NES_mat, PES_mat]
    return result
