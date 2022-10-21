'''
Module for selecting individuals with large breeding values
'''

import sys
from copy import deepcopy
from random import shuffle
from operator import itemgetter
import numpy as np
from aquabreeding import blup as bl


def select_value(pro_inf, phe_inf, target):
    '''
    Set values, by which progenies are selected

    Args:
        pro_inf (PopulationInfo class): Progeny population
        phe_inf (PhenotypeInfo class): Phenotype information
        target (str): 'bv' (breeding value), 'phenotype' (phenotype), or
                      'random' (random) selection

    Returns:
        ndarray: breeding value, phenotype or random number
    '''
    if target == 'bv':
        return phe_inf.hat_bv
    if target == 'phenotype':
        return phe_inf.pheno_v
    if target == 'random':
        n_founder = pro_inf.n_popf + pro_inf.n_popm
        return np.random.rand(n_founder)
    sys.exit('target should be \'bv\', \'phenotype\', or \'random\'')
# select_value


def family_select_value(pro_inf, select_val, fam_d, tmp_bv, f_m):
    '''
    Get index, family, and, breeding value

    Args:
        pro_inf (PopulationInfo class): Progeny population
        select_val (ndarray): Breeding value or others
        fam_d (dict): Dictionary for family
        tmp_bv (list): Index, family, value, 0 (female) or 1 (male)
        f_m (int): 0 for female, 1 for male
    '''
    if f_m == 0:  # for female
        p_pop = pro_inf.pop_f
        n_pop = pro_inf.n_popf
        sv_i = 0
    elif f_m == 1:  # for male
        p_pop = pro_inf.pop_m
        n_pop = pro_inf.n_popm
        sv_i = pro_inf.n_popf
    else:
        sys.exit('0 for female or 1 for male')
    for i in range(n_pop):
        if p_pop[i].pat_id < p_pop[i].mat_id:
            fam_n = f'f{p_pop[i].pat_id}x{p_pop[i].mat_id}'
        else:
            fam_n = f'f{p_pop[i].mat_id}x{p_pop[i].pat_id}'
        tmp_bv.append([i, fam_n, select_val[i+sv_i], f_m])
        fam_d.setdefault(fam_n, 0)
# family_select_value


def get_index(par_inf, pro_inf, tmp_bv, fam_d, top_prop):
    '''
    Get Index or selected progenies

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        tmp_bv (list): index, family, value, 0 (female) or 1 (male).
                       Sorted by value.
        fam_d (dict): Dictionary of family (all vales are zero)
        top_prop (float): Select progenies with top X% values

    Returns:
        Index of selected individuals


        - ndarray: Index of females
        - ndarray: Index of males
    '''
    f_val = np.full(par_inf.n_popf, -7, dtype=np.int64)
    m_val = np.full(par_inf.n_popm, -7, dtype=np.int64)
    n_fval = n_mval = 0
    n_tolerant = 0  # how many times progenies are selected in a family
    n_p2 = int(top_prop * (pro_inf.n_popf + pro_inf.n_popm))
    if n_p2 < (par_inf.n_popf + par_inf.n_popm):
        sys.exit('top_prop is too small')
    while n_fval != par_inf.n_popf or n_mval != par_inf.n_popm:
        for i in range(n_p2):
            if tmp_bv[i][0] == -7:
                continue
            if fam_d[tmp_bv[i][1]] != n_tolerant:
                continue
            if tmp_bv[i][3] == 0:  # female
                if n_fval == par_inf.n_popf:
                    continue
                f_val[n_fval] = tmp_bv[i][0]
                n_fval += 1
            else:  # male
                if n_mval == par_inf.n_popm:
                    continue
                m_val[n_mval] = tmp_bv[i][0]
                n_mval += 1
            tmp_bv[i][0] = -7
            fam_d[tmp_bv[i][1]] += 1
        n_tolerant += 1
        if n_tolerant >= pro_inf.n_popf + pro_inf.n_popm:
            sys.exit('top_prop or n_family is too small')
    return f_val, m_val
# get_index


def within_family_selection(par_inf, pro_inf, select_val, top_prop):
    '''
    Get index of progenies by within-family selection

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        select_val (ndarray): Selection based on these values
        top_prop (float): Select progenies with top X% of values

    Returns:
        Index of selected individuals

        - ndarray: Index of females
        - ndarray: Index of males
    '''
    fam_d = {}  # pair of mother and father
    tmp_bv = []  # list of family, gender, and bv
    family_select_value(pro_inf, select_val, fam_d, tmp_bv, 0)  # for female
    family_select_value(pro_inf, select_val, fam_d, tmp_bv, 1)  # for fale
    # sort by breeding value (or others)
    # [0: ID, 1: family, 2: breeding value, 3: female or male]
    tmp_bv.sort(reverse=True, key=itemgetter(2))
    f_val, m_val = get_index(par_inf, pro_inf, tmp_bv, fam_d, top_prop)
    return f_val, m_val
# within_family_selection


def family_mean_value(tmp_bv, fam_d, n_family):
    '''
    Get mean values of each family

    Args:
        tmp_bv (list): Index, family, value, 0 (female) or 1 (male)
        fam_d (dict): Dictionary for family
    '''
    # check if n_family is too large
    if len(fam_d) < n_family:
        n_family = len(fam_d)
    tmp_d = {}
    tmp_fam = []
    # initialize
    for i in fam_d.keys():
        tmp_d[i] = []
    # get values for each family
    for i in tmp_bv:
        tmp_d[i[1]].append(i[2])
    # list[family, mean value]
    for k, v_ls in tmp_d.items():
        tmp_fam.append([k, np.mean(v_ls)])
    tmp_fam.sort(reverse=True, key=itemgetter(1))
    # set of used family
    used_fam = []
    for i in range(n_family):
        used_fam.append(tmp_fam[i][0])
    used_fam = set(used_fam)
    # mask unused family
    for i in tmp_bv:
        if not i[1] in used_fam:
            i[0] = -7
# family_mean_value


def family_selection(par_inf, pro_inf, select_val, n_family):
    '''
    Get index of progenies by family selection

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        select_val (ndarray): Selection based on these values
        n_family (int): Number of families to be selected

    Returns:
        Index of selected individual

        - ndarray: Index of females
        - ndarray: Index of males
    '''
    if n_family == -1:
        sys.exit('n_family should be set to use family selection')
    fam_d = {}  # pair of mother and father
    tmp_bv = []  # list of family, gender, and bv
    family_select_value(pro_inf, select_val, fam_d, tmp_bv, 0)  # for female
    family_select_value(pro_inf, select_val, fam_d, tmp_bv, 1)  # for fale
    # sort by breeding value (or others)
    # [0: ID, 1: family, 2: breeding value, 3: female or male]
    tmp_bv.sort(reverse=True, key=itemgetter(2))
    # mask unused families
    family_mean_value(tmp_bv, fam_d, n_family)
    # get index
    f_val, m_val = get_index(par_inf, pro_inf, tmp_bv, fam_d, 1.0)
    return f_val, m_val
# family_selection


def mass_selection(par_inf, pro_inf, select_val):
    '''
    Get index of progenies with the largest select values

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        select_val (ndarray): Selection based on these values

    Returns:
        Index of selected individuals

        - ndarray: index of females
        - ndarray: Index of males
    '''
    f_val = np.argsort(select_val[:pro_inf.n_popf])[::-1][:par_inf.n_popf]
    m_val = np.argsort(select_val[pro_inf.n_popf:])[::-1][:par_inf.n_popm]
    return f_val, m_val
# mass_selection


def nextgen_parents(par_inf, pro_inf, f_val, m_val, cross_inf):
    '''
    Set next-generation parents as inbreeding is minimized

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        f_val (ndarray): Index of selected females
        m_val (ndarray): Index of selected males
        cross_inf (ndarray): Female index, male index, no. female
                             progeny, no. male progeny
    '''
    # numerator relationship matrix of selected ones
    n_mat = bl.nrm_selected_ones(par_inf, pro_inf, f_val, m_val)
    # random search for minimum mean and max IBD
    res_ibd = []
    for i in range(1000):
        f_1 = list(range(par_inf.n_popf))
        m_1 = list(range(par_inf.n_popm))
        shuffle(f_1)
        shuffle(m_1)
        ibd_tmp = []
        for k_1, k_2, k_3, k_4 in cross_inf:
            ibd_tmp.append(n_mat[f_1[k_1]][m_1[k_2]+par_inf.n_popf])
        mean_ibd = np.mean(ibd_tmp)
        max_ibd = max(ibd_tmp)
        res_ibd.append([f_1, m_1, mean_ibd, max_ibd])
    # get pairs with minimum inbreeding
    res_ibd.sort(key=itemgetter(2, 3))
    f_2 = res_ibd[0][0]
    m_2 = res_ibd[0][1]
    # female
    for i in range(par_inf.n_popf):
        pro_i = f_val[f_2[i]]
        par_inf.pop_f[i] = deepcopy(pro_inf.pop_f[pro_i])
    # male
    for i in range(par_inf.n_popm):
        pro_i = m_val[m_2[i]]
        par_inf.pop_m[i] = deepcopy(pro_inf.pop_m[pro_i])
    # add new parents' ID
    par_inf.new_founder_id()
# nextgen_parent


def select_parent(par_inf, pro_inf, phe_inf, cross_inf, target, method,
                  top_prop, n_family):
    '''
    Select parent

    Args:
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        phe_inf (PhenotypeInfo class): Phenotype information
        cross_inf (ndarray): Female index, male index, no. female progeny,
                             no. male progeny
        target (str): 'bv' (breeding value), 'phenotype' (phenotype), or
                       'random' (random) selection
        method (str): 'mass' (mass selection), 'within-family' (within-
                      family selection), or 'family' (family selection) 
        top_prop (float): Individual with top X% breeding values are used
                          in within-family selection
        n_family (int): Number of families to be selected in family selection
    '''
    # breeding value, phenotype, or random
    select_val = select_value(pro_inf, phe_inf, target)
    # get ID of large select values
    # mass selection
    if method == 'mass':
        f_val, m_val = mass_selection(par_inf, pro_inf, select_val)
    # within-family selection
    elif method == 'within-family':
        f_val, m_val = within_family_selection(par_inf, pro_inf, select_val,
                                               top_prop)
    # family selection
    elif method == 'family':
        f_val, m_val = family_selection(par_inf, pro_inf, select_val, n_family)
    else:
        sys.exit('method should be \'mass\', or \'within-family\'')
    # mating design as inbreeding is minimized
    nextgen_parents(par_inf, pro_inf, f_val, m_val, cross_inf)
# select_parent


def main():
    '''
    main
    '''
    print('Module for selection')
# main


if __name__ == '__main__':
    main()
