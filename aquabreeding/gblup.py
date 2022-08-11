'''
module for GBLUP
'''

import sys
import numpy as np
import scipy
from numba import jit
from aquabreeding import _nrm as nrm


def nrm_cpp(par_dict, progeny_pop, founder_size):
    '''
    version 3
    calculating numerator relationship matrix
    parental relationship of all founders and current progeny are required
    args: dictionary of parents of founders in all generations
          current progeny populatin info
          size of founder population
    return np.ndarray (n x n)
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


def gebv_calculation(y_vec, g_mat):
    '''
    version 4
    Calculating breeding values and variance components
    by a restricted maximum likelihood method
    The script is converted from the R package, rrBLUP, by
    Endelman (2011) Plant Genome 4, 250-255.
    The original method is described in
    Kang et al. (2008) Genetics 178, 1709-1723.
    '''
    def restricted_ml(lam, n_p, theta_vec, omega_sq):
        return (n_p*np.log(np.sum(omega_sq/(theta_vec + lam))) +
                np.sum(np.log(theta_vec + lam)))
    n_sample = y_vec.shape[0]
    # design matrix of X and Z
    x_mat = np.ones((n_sample, 1))
    z_mat = np.identity(n_sample)
    # the number of fixed effects
    n_p = x_mat.shape[1]
    # Kang et al. (2008)
    xt_mat = x_mat.T
    xtx_mat = xt_mat @ x_mat
    xtx_inv = np.linalg.inv(xtx_mat)
    s_mat = np.identity(n_sample) - x_mat @ xtx_inv @ xt_mat
    off_set = np.sqrt(n_sample)
    hb_mat = z_mat @ g_mat @ z_mat.T + off_set*np.identity(n_sample)
    hb_val, hb_vec = np.linalg.eigh(hb_mat)
    phi_vec = hb_val.ravel() - off_set
    s_hb_s = s_mat @ hb_mat @ s_mat
    shbs_val, shbs_vec = np.linalg.eigh(s_hb_s)
    theta_vec = shbs_val[n_p:] - off_set
    q_mat = shbs_vec[:, n_p:]
    omega_sq = ((q_mat.T @ y_vec)**2).ravel()
    set_arg = (n_sample-n_p, theta_vec, omega_sq)
    opt_res = scipy.optimize.minimize_scalar(restricted_ml,
                                             args=set_arg,
                                             bounds=(1e-9, 1e9),
                                             method='bounded')
    if not opt_res.success:
        sys.exit('scipy.optimize does not work')
    opt_para = opt_res.x
    d_f = n_sample - n_p
    # estimated variance components
    hat_vg = np.sum(omega_sq/(theta_vec + opt_para))/d_f
    hat_ve = opt_para*hat_vg
    h_inv = hb_vec @ (hb_vec/(phi_vec+opt_para)).T
    xt_hinv = xt_mat @ h_inv
    w_mat = xt_hinv @ x_mat
    # estimated fixed effects
    hat_beta = np.linalg.solve(w_mat, xt_hinv @ y_vec)
    kz_t = g_mat @ z_mat.T
    kzt_hinv = kz_t @ h_inv
    # predicted breeding value
    hat_bv = kzt_hinv @ (y_vec - x_mat @ hat_beta)
    return hat_beta, hat_bv, hat_vg, hat_ve
    # gebv_calculation


@jit
def nrm_jit2(n_mat, j_stt, j_max, vec_s, vec_d):
    '''
    calculating numerator relationship matrix
    with numba jit
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


def nrm_jit(par_dict, pro_pop, n_founder):
    '''
    version 4
    calculating numerator relationship matrix
    with numba jit
    '''
    par_key = list(par_dict.keys())
    par_key.sort()
    n_par = len(par_key)
    dict_s = []
    dict_d = []
    for i in range(n_par):
        dict_s.append(par_dict[i]['pat'])
        dict_d.append(par_dict[i]['mat'])
        n_progeny = len(pro_pop)
    for i in range(n_progeny):
        dict_s.append(pro_pop[i].pat_id)
        dict_d.append(pro_pop[i].mat_id)
    n_mat = np.identity(n_par+n_progeny)
    dict_s = np.array(dict_s)
    dict_d = np.array(dict_d)
    nrm_jit2(n_mat, n_founder, n_par+n_progeny, dict_s, dict_d)
    return n_mat[n_par:, n_par:]
    # nrm_jit


def main():
    '''
    main
    '''
    print('The module for GBLUP')
    # main


if __name__ == '__main__':
    main()
