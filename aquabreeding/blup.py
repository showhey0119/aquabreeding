'''
Module for BLUP
'''

import sys
from scipy.optimize import minimize_scalar
import numpy as np
from numba import jit
from aquabreeding import _nrm as nrm


def nrm_cpp(par_dict, progeny_pop, founder_size):
    '''
    Calculate numerator relationship matrix with C++

    Args:
        par_dict (dict): Dictionary of parents in founders
        progeny_pop (PopulationInfo class): Progeny population
        founder_size (int): Founder population size

    Returns:
        (ndarray: Numerator relationship matrix
    '''
    par_key = list(par_dict.keys())
    par_key.sort()
    n_par = len(par_key)
    dict_s = []
    dict_d = []
    for i in range(n_par):
        dict_s.append(par_dict[i]['pat'])
        dict_d.append(par_dict[i]['mat'])
    progeny_size = len(progeny_pop)
    for i in range(progeny_size):
        dict_s.append(progeny_pop[i].pat_id)
        dict_d.append(progeny_pop[i].mat_id)
    n_mat = np.identity(n_par+progeny_size)
    nrm.nrm(n_mat, founder_size, n_par+progeny_size, dict_s, dict_d)
    return n_mat[n_par:, n_par:]
# nrm


def bv_estimation(phe_inf, g_mat):
    '''
    Estimate breeding values and variance components
    Version 5:
    Estimate breeding values and variance components
    by a restricted maximum likelihood method.
    The script is converted from the R package, rrBLUP, by
    Endelman (2011) Plant Genome 4, 250-255.
    The original method is described in
    Kang et al. (2008) Genetics 178, 1709-1723.
    Args:
        phe_inf (PhenotypeInfo class): Phenotype information
        g_mat (ndarray): Genomic/numerator relationship matrix
    Note:
        Fixed effect is not implement yet.
    '''
    def restricted_ml(lam, n_p, theta_vec, omega_sq):
        return (n_p*np.log(np.sum(omega_sq/(theta_vec + lam))) +
                np.sum(np.log(theta_vec + lam)))
    y_vec = phe_inf.pheno_v
    n_sample = y_vec.shape[0]
    # design matrix of X and Z
    x_mat = np.ones((n_sample, 1), dtype=np.float64)
    u_1, s_1, vh_1 = np.linalg.svd(x_mat)
    r_i = int(max(np.where(s_1 > 1e-8)))
    if r_i == 0:
        x_mat = u_1[:,r_i].reshape(n_sample, 1)
    else:
        sys.exit('Something wrong in bv_estimation')
    z_mat = np.identity(n_sample, dtype=np.float64)
    # the number of fixed effects
    n_p = x_mat.shape[1]
    # Convert rrBLUP/R into python
    xt_mat = x_mat.T
    xtx_mat = xt_mat @ x_mat
    xtx_inv = np.linalg.inv(xtx_mat)
    s_mat = np.identity(n_sample, dtype=np.float64) - x_mat @ xtx_inv @ xt_mat
    off_set = np.sqrt(n_sample)
    hb_mat = z_mat @ g_mat @ z_mat.T + off_set * np.identity(n_sample)
    hb_val, hb_vec = np.linalg.eigh(hb_mat)
    phi_vec = hb_val.ravel() - off_set
    s_hb_s = s_mat @ hb_mat @ s_mat
    shbs_val, shbs_vec = np.linalg.eigh(s_hb_s)
    theta_vec = shbs_val[n_p:] - off_set
    q_mat = shbs_vec[:, n_p:]
    omega_sq = ((q_mat.T @ y_vec)**2).ravel()
    set_arg = (n_sample-n_p, theta_vec, omega_sq)
    opt_res = minimize_scalar(restricted_ml, args=set_arg, bounds=(1e-9, 1e9),
                              method='bounded')
    if not opt_res.success:
        sys.exit('scipy.optimize does not work')
    opt_para = opt_res.x
    d_f = n_sample - n_p
    # estimate variance components
    phe_inf.hat_vg = np.sum(omega_sq/(theta_vec + opt_para))/d_f
    phe_inf.hat_ve = opt_para*phe_inf.hat_vg
    h_inv = hb_vec @ (hb_vec/(phi_vec+opt_para)).T
    xt_hinv = xt_mat @ h_inv
    w_mat = xt_hinv @ x_mat
    # estimate fixed effects
    phe_inf.hat_beta = np.linalg.solve(w_mat, xt_hinv @ y_vec)    
    kz_t = g_mat @ z_mat.T
    kzt_hinv = kz_t @ h_inv
    # estimate breeding value
    phe_inf.hat_bv = kzt_hinv @ (y_vec - x_mat @ phe_inf.hat_beta)
    phe_inf.hat_beta = x_mat.mean(axis=0) @ phe_inf.hat_beta
# gebv_calculation


@jit(cache=True)
def nrm_jit2(n_mat, j_stt, j_max, vec_s, vec_d):
    '''
    Calculating numerator relationship matrix with numba jit

    Args:
        n_mat (ndarray): Numerator relationship matrix
        j_stt (int): Founder size
        j_max (int): The total number of founder+progeny
        vec_s (ndarray): Father (sire) ID
        vec_d (ndarray): Mother (dam) ID
    '''
    for j in range(j_stt, j_max):
        pars = vec_s[j]
        pard = vec_d[j]
        for i in range(j):
            ais = 0.0
            aid = 0.0
            if pars != -1:
                ais = n_mat[i][pars]
            if pard != -1:
                aid = n_mat[i][pard]
            n_mat[i][j] = 0.5*(ais+aid)
            n_mat[j][i] = n_mat[i][j]
        if pars != -1 and pard != -1:
            n_mat[j][j] = 1.0 + 0.5 * n_mat[pars][pard]
        else:
            n_mat[j][j] = 1.0
# nrm_jit2


def nrm_jit(par_inf, pro_inf):
    '''
    Calculate numerator relationship matrix

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population

    Returns:
        ndarray: Numerator relationship matrix
    '''
    # parents of the founder population
    par_key = list(par_inf.d_par.keys())
    par_key.sort()
    n_founder = par_inf.n_popf + par_inf.n_popm
    # num of parents throughout all generations
    n_par = len(par_key)
    # list of father and mother
    dict_s = []
    dict_d = []
    for i in range(n_par):
        dict_s.append(par_inf.d_par[i]['pat'])
        dict_d.append(par_inf.d_par[i]['mat'])
    # parents of the progeny population
    n_progeny = pro_inf.n_popf + pro_inf.n_popm
    for i in range(pro_inf.n_popf):
        dict_s.append(pro_inf.pop_f[i].pat_id)
        dict_d.append(pro_inf.pop_f[i].mat_id)
    for i in range(pro_inf.n_popm):
        dict_s.append(pro_inf.pop_m[i].pat_id)
        dict_d.append(pro_inf.pop_m[i].mat_id)
    dict_s = np.array(dict_s, dtype=np.int64)
    dict_d = np.array(dict_d, dtype=np.int64)
    # calculate numerator relationship matrix
    n_mat = np.identity(n_par+n_progeny, dtype=np.float64)
    nrm_jit2(n_mat, n_founder, n_par+n_progeny, dict_s, dict_d)
    return n_mat[n_par:, n_par:]
# nrm_jit


def nrm_selected_ones(par_inf, pro_inf, f_val, m_val):
    '''
    Calculate numerator relationship matrix of selected ones

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        f_val (ndarray): Female index of selected progenies
        m_val (ndarray): Male index of selected progenies

    Returns:
        ndarray: Numerator relationship matrix
    '''
    # parents of the founder population
    par_key = list(par_inf.d_par.keys())
    par_key.sort()
    n_founder = par_inf.n_popf + par_inf.n_popm
    # num of parents throughout all generations
    n_par = len(par_key)
    # list of father and mother
    dict_s = []
    dict_d = []
    for i in range(n_par):
        dict_s.append(par_inf.d_par[i]['pat'])
        dict_d.append(par_inf.d_par[i]['mat'])
    # parents of selected progeny
    n_progeny = f_val.shape[0] + m_val.shape[0]
    for i in f_val:
        dict_s.append(pro_inf.pop_f[i].pat_id)
        dict_d.append(pro_inf.pop_f[i].mat_id)
    for i in m_val:
        dict_s.append(pro_inf.pop_m[i].pat_id)
        dict_d.append(pro_inf.pop_m[i].mat_id)
    dict_s = np.array(dict_s)
    dict_d = np.array(dict_d)
    # calculate numerator relationship matrix
    n_mat = np.identity(n_par+n_progeny)
    nrm_jit2(n_mat, n_founder, n_par+n_progeny, dict_s, dict_d)
    return n_mat[n_par:, n_par:]
# nrm_selected_ones


def ablup(phe_inf, par_inf, pro_inf):
    '''
    BLUP with numerator relationship matrix (ABLUP)

    Args:
        phe_inf (PhenotypeInfo class): Phenotype information
        par_inf (PopulationInfo class): Founder poplation
        pro_inf (PopulationInfo class): Progeny population
    '''
    # numerator relationship matrix
    a_mat = nrm_jit(par_inf, pro_inf)
    # estimate breeding value and variance component
    bv_estimation(phe_inf, a_mat)
# ablup


@jit(cache=True)
def convert_gmatrix(gen_array, p_mat, gmat_denom):
    '''
    Convert genotype matrix into G matrix for GBLUP

    Args:
        gen_array (ndarray): Genotype matrix
        p_mat (ndarray): Allele frequency matrix (P matrix)
        gmat_denom (float): Normalize G matrix by this

    Returns:
        ndarray: G matrix
    '''
    n_cols = gen_array.shape[1]
    # calculate denominator for normalizing the matrix
    for i in range(n_cols):
        # i th locus
        gen_array[:, i] -= 2.0 * p_mat[i]
    g_mat = gen_array @ gen_array.T / gmat_denom
    return g_mat
# convert_gmatrix


def gblup(phe_inf, gblup_inf):
    '''
    BLUP with genomic relationship matrix (GBLUP)

    Args:
        phe_inf (PhenotypeInfo class): Phenotype information
        gblup_inf (SNPInfo class): SNP information for GBLUP
    '''
    # calculate G matrix from genotype matrix
    g_mat = convert_gmatrix(gblup_inf.gen_mat, gblup_inf.p_mat,
                            gblup_inf.gmat_denom)
    # estimate breeding value and variance component
    bv_estimation(phe_inf, g_mat)
# gblup


def main():
    '''
    main
    '''
    print('Module for BLUP')
# main


if __name__ == '__main__':
    main()
