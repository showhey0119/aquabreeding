'''
A module for BLUP
'''

import sys
from scipy.optimize import minimize_scalar
import numpy as np
from aquabreeding import aquacpp as cpp


def nrm_cpp(par_inf, pro_inf):
    '''
    Calculate numerator relationship matrix

    Args:
        par_inf (PopulationInfo): Founder population
        pro_inf (PopulationInfo): Progeny population
    '''
    par_dict = par_inf.d_par
    par_key = list(par_dict.keys())
    par_key.sort()
    n_par = len(par_key)
    dict_s = []
    dict_d = []
    for i in range(n_par):
        dict_s.append(par_dict[i]['pat'])
        dict_d.append(par_dict[i]['mat'])
    # female
    for individual in pro_inf.pop_f:
        dict_s.append(individual.pat_id)
        dict_d.append(individual.mat_id)
    # male
    for individual in pro_inf.pop_m:
        dict_s.append(individual.pat_id)
        dict_d.append(individual.mat_id)
    n_founder = par_inf.n_founder
    n_progeny = pro_inf.n_f + pro_inf.n_m
    n_mat = np.identity(n_par + n_progeny, dtype=np.float64)
    cpp.nrm(n_mat, n_founder, n_par + n_progeny, dict_s, dict_d)
    pro_inf.a_mat = n_mat[n_par:, n_par:]
# nrm_cpp


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
        phe_inf (PhenotypeInfo): Phenotype information
        g_mat (numpy.ndarray): Genomic/numerator relationship matrix
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
        x_mat = u_1[:, r_i].reshape(n_sample, 1)
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
# bv_estimation


def convert_gmatrix(pro_inf, gen_array):
    '''
    Convert genotype matrix into G matrix for GBLUP

    Args:
        pro_inf (PopulationInfo): Progeny population
        gen_array (numpy.ndarray): Genotype matrix

    Returns:
        numpy.ndarray: Genomic relationship matrix
    '''
    n_cols = gen_array.shape[1]
    # allele frequency
    freq_q = np.ones(n_cols)
    for i in range(n_cols):
        freq_q[i] = np.mean(gen_array[:, i])/2.0
    # Remove monomorphic
    poly_snp = freq_q * (1.0 - freq_q) > 1e-10
    freq_q = freq_q[poly_snp]
    # Genomic relationship matrix
    var_a = 2 * np.mean(freq_q * (1.0 - freq_q))
    w_mat = gen_array[:, poly_snp] - 2.0 * freq_q
    n_poly = len(freq_q)
    g_mat = w_mat @ w_mat.T / var_a / n_poly
    pro_inf.g_mat = g_mat
# convert_gmatrix


def main():
    '''
    main
    '''
    print('A module for BLUP')
# main


if __name__ == '__main__':
    main()
