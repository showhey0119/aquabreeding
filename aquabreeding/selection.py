'''
A module for selection
'''

import sys
from math import isclose
from random import choice
from operator import itemgetter
import numpy as np
from aquabreeding import aquacpp as cpp


def select_value(target, phe_inf):
    '''
    Set values, by which progenies are selected

    Args:
        target (str): 'bv' (breeding value), 'phenotype' (phenotype), or
                      'random' (random) selection
        phe_inf (PhenotypeInfo): Phenotype information

    Returns:
        numpy.ndarray: breeding value, phenotype or random number
    '''
    if target == 'bv':
        return phe_inf.hat_bv
    if target == 'phenotype':
        return phe_inf.pheno_v
    if target == 'random':
        n_pro = np.shape(phe_inf.pheno_v)[0]
        return np.random.rand(n_pro)
    sys.exit('Invalid in select_value')
# select_value


def within_family_selection(summary_pro, select_size, fam_d, top_prop):
    '''
    Get index of progenies by within-family selection

    Args:
        summary_pro (list): Index, family, select value, female or male
        select_size (tuple): Nos. selected parents
        fam_d (dict): family dict
        top_prop (float): Select progenies from top X% values

    Returns:
        Index of selected individuals

        - numpy.ndarray: Index of females
        - numpy.ndarray: Index of males
    '''
    f_index = np.full(select_size[0], -7, dtype=np.int32)
    m_index = np.full(select_size[1], -7, dtype=np.int32)
    n_f = n_m = 0
    n_tolerant = 0  # how many times progenies are selected in a family
    n_p2 = int(top_prop * len(summary_pro))
    if n_p2 < (select_size[0] + select_size[1]):
        sys.exit('top_prop is too small')
    while n_f != select_size[0] or n_m != select_size[1]:
        for i in range(n_p2):
            tmp_ls = summary_pro[i]
            # already selected or masked
            if tmp_ls[0] == -7:
                continue
            # Enough no. individuals are selected in the family
            if fam_d[tmp_ls[1]] != n_tolerant:
                continue
            # female
            if tmp_ls[3] == 0 and n_f < select_size[0]:
                f_index[n_f] = tmp_ls[0]
                n_f += 1
            # male
            elif tmp_ls[3] == 1 and n_m < select_size[1]:
                m_index[n_m] = tmp_ls[0]
                n_m += 1
            else:
                continue
            tmp_ls[0] = -7
            fam_d[tmp_ls[1]] += 1
        n_tolerant += 1
        if n_tolerant >= select_size[0] + select_size[1]:
            sys.exit('top_prop or n_family is too small')
    return f_index, m_index
# within_family_selection


def mass_selection(summary_pro, select_size):
    '''
    Get index of parents by mass selection

    Args:
        summary_pro (list): Index, family, select value, female or male
        select_size (tuple): Nos. selected parents
    '''
    f_index = np.full(select_size[0], -7, dtype=np.int32)
    m_index = np.full(select_size[1], -7, dtype=np.int32)
    n_f = n_m = 0
    for tmp_ls in summary_pro:
        # female
        if tmp_ls[3] == 0 and n_f < select_size[0]:
            f_index[n_f] = tmp_ls[0]
            n_f += 1
        # male
        elif tmp_ls[3] == 1 and n_m < select_size[1]:
            m_index[n_m] = tmp_ls[0]
            n_m += 1
        # end selection
        if n_f == select_size[0] and n_m == select_size[1]:
            return f_index, m_index
    return -7, -7
# mass_selection


def merge_pro_info(pop_x, select_val, fam_d, summary_pro, i_rows, f_m):
    '''
    Merge progeny info

    Args:
        pop_x (list): List of IndividualInfo class
        select_val (numpy.ndarray): Breeding value or others
        fam_d (dict): Dictionary for family
        summary_pro (list): Index, family, value, 0 (female) or 1 (male)
        i_rows (int): Start index
        f_m (int): 0 for female, 1 for male
    '''
    for i, individual in enumerate(pop_x):
        # family name
        fam_n = f'f{individual.mat_id}x{individual.pat_id}'
        fam_d.setdefault(fam_n, 0)
        # index, family, select value, female or male
        summary_pro.append([i, fam_n, select_val[i+i_rows], f_m])
    summary_pro.sort(reverse=True, key=itemgetter(2))
# merge_pro_info


def mask_family(summary_pro, fam_d, n_family):
    '''
    Mask unused family

    Args:
        summary_pro (list): Index, family, select value, female or male
        fam_d (dict): Family dict
        n_family (int): No. families to be selected
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
    for i in summary_pro:
        tmp_d[i[1]].append(i[2])
    # list[family, mean value]
    for k, v_ls in tmp_d.items():
        tmp_fam.append([k, np.mean(v_ls)])
    tmp_fam.sort(reverse=True, key=itemgetter(1))
    # set of used family
    used_fam = set()
    for i in range(n_family):
        used_fam.add(tmp_fam[i][0])
    # mask unused family
    for i in summary_pro:
        if not i[1] in used_fam:
            i[0] = -7
# mask_family


def family_selection(summary_pro, select_size, fam_d, top_prop, n_family):
    '''
    Get index of progenies by family selection

    Args:
        summary_pro (list): Index, family select value, female or male
        select_size (tuple): Nos. selected parents
        fam_d (dict): family dict
        top_prop (float): Select progenies with top X% values
        n_family (int): No. families to be selected

    Returns:
        Index of selected individuals

        - numpy.ndarray: Index of females
        - numpy.ndarray: Index of males
    '''
    if n_family == -1:
        sys.exit('n_family should be set to use family selection')
    # mask unused families
    mask_family(summary_pro, fam_d, n_family)
    # get index
    f_index, m_index = within_family_selection(summary_pro, select_size, fam_d,
                                               top_prop)
    return f_index, m_index
# family_selection


def fvalue_selection(summary_pro, select_size, top_prop, g_mat, max_r, i_rows):
    '''
    Get Index or selected progenies based on G matrix

    Args:
        summary_pro (list): Index, family, select value, female or male
        select_size (tuple): Nos. parents to be selected
        top_prop (float): Select progenies from top X% values
        g_mat (numpy.ndarray): Numerator/genomic relationship matrix
        max_r (float): Max F value
        i_rows (int): To get male F value from A/G matrix

    Returns:
        Index of selected individuals

        - numpy.ndarray: Index of females
        - numpy.ndarray: Index of males
    '''
    f_index = np.full(select_size[0], -7, dtype=np.int32)
    m_index = np.full(select_size[1], -7, dtype=np.int32)
    n_f = n_m = 0
    n_bv = int(top_prop * len(summary_pro))
    for i in range(n_bv):
        tmp_ls = summary_pro[i]
        if tmp_ls[0] == -7:
            continue
        # G matrix index
        i_1 = tmp_ls[0]
        # famale
        if tmp_ls[3] == 0:
            if n_f == select_size[0]:
                continue
            f_index[n_f] = tmp_ls[0]
            n_f += 1
        else:  # male
            if n_m == select_size[1]:
                continue
            m_index[n_m] = tmp_ls[0]
            n_m += 1
            i_1 += i_rows
        tmp_ls[0] = -7
        # Mask if F value is larger than threshold
        for j in range(i+1, n_bv):
            tmp_ls2 = summary_pro[j]
            # G matrix index
            i_2 = tmp_ls2[0]
            if tmp_ls2[3] == 1:
                i_2 += i_rows
            # If two individuals are closely related
            if g_mat[i_1][i_2] > max_r:
                tmp_ls2[0] = -7
    if n_f != select_size[0] or n_m != select_size[1]:
        print(f'Not enough {select_size}, {n_f} {n_m}, {max_r:.1f}')
        return f_index, m_index, 1
        sys.exit(f'Not enough {select_size}, {n_f} {n_m}, {max_r:.2f}')
    return f_index, m_index, 0
# fvalue_selection


def get_furthest_one(g_mat, x_index, y_ls, i_rows, f_m):
    '''
    Get the most distantly related one

    Args:
        g_mat (numpy.ndarray): Numerator/Genomic relationship matrix
        x_index (int): female/male index
        y_ls (list): remained male/female index list
        i_rows (int): To get index of male in G matrix
        f_m (int): 0: x_index is female, 1: male

    Returns:
        int: index of distantly related one
    '''
    tmp_ls = []
    # x_index for female
    if f_m == 0:
        for y_1 in y_ls:
            tmp_ls.append([y_1, g_mat[x_index][y_1+i_rows]])
    else:
        for y_1 in y_ls:
            tmp_ls.append([y_1, g_mat[x_index+i_rows][y_1]])
    # minimum inbreeding coefficient
    min_ic = min(i[1] for i in tmp_ls)
    tmp_ls2 = [i[0] for i in tmp_ls if isclose(i[1], min_ic)]
    pick_one = choice(tmp_ls2)
    y_ls.remove(pick_one)
    return pick_one
# get_furthest_one


def arrange_parents(g_mat, cross_info, select_size, f_index, m_index,
                    i_rows):
    '''
    Arrange selected parents based on mating design

    Args:
        g_mat (numpy.ndarray): Numerator/genomic relationship  matrix
        cross_inf (numpy.ndaray): Index pairs of female and male parents
        select_size (tuple): No. parents to be selected
        f_index (numpy.ndarray): Index of selected females
        m_index (numpy.ndarray): Index of selected males
        i_rows (int): To get index of males in G matrix

    Returns:
        List of index

        - list: Female index
        - list: Male index
    '''
    n_f, n_m = select_size
    # serach good arrangement
    res_ibd = []
    for _ in range(100):
        # key: index in parental pop list, value: selected progeny index
        f_dict = {i: -1 for i in range(n_f)}
        m_dict = {i: -1 for i in range(n_m)}
        # remained progeny index, value = -1
        f_ls = list(f_index)
        m_ls = list(m_index)
        # arrange progenies as inbreeding is minimized
        for c_1, c_2 in cross_info:
            # no info
            if f_dict[c_1] == -1 and m_dict[c_2] == -1:
                # get random female
                f_dict[c_1] = choice(f_ls)
                f_ls.remove(f_dict[c_1])
                # get male with no shared ancestor
                m_dict[c_2] = get_furthest_one(g_mat, f_dict[c_1], m_ls,
                                               i_rows, 0)
            elif f_dict[c_1] != -1 and m_dict[c_2] == -1:
                m_dict[c_2] = get_furthest_one(g_mat, f_dict[c_1], m_ls,
                                               i_rows, 0)
            elif f_dict[c_1] == -1 and m_dict[c_2] != -1:
                f_dict[c_1] = get_furthest_one(g_mat, m_dict[c_2], f_ls,
                                               i_rows, 1)
        f_1 = [f_dict[i] for i in range(n_f)]
        m_1 = [m_dict[i] for i in range(n_m)]
        # check inbreeding
        ibd_tmp = []
        for c_1, c_2 in cross_info:
            ibd_tmp.append(g_mat[f_dict[c_1]][m_dict[c_2]+i_rows])
        mean_ibd = np.mean(ibd_tmp)
        max_ibd = max(ibd_tmp)
        res_ibd.append([f_1, m_1, mean_ibd, max_ibd])
    # choose the best
    res_ibd.sort(key=itemgetter(2, 3))
    f_2 = res_ibd[0][0]
    m_2 = res_ibd[0][1]
    return f_2, m_2
# arrange_parents


def copy_individual(ls_founder, ls_progeny, ls_index):
    '''
    Copy individual information

    Args:
        ls_founder (list): List of IndividualInfo class in founder
        ls_progeny (list): List of IndividualInfo class in progeny
        ls_index (list): List of index of selected progeny
    '''
    for founder, x_index in zip(ls_founder, ls_index):
        par_chrom_ls = founder.chrom_ls
        pro_chrom_ls = ls_progeny[x_index].chrom_ls
        for par_chrom, pro_chrom in zip(par_chrom_ls, pro_chrom_ls):
            if par_chrom.position is None:
                par_chrom.position = pro_chrom.position.copy()
            else:
                cpp.copy_1D(pro_chrom.position, par_chrom.position)
            if par_chrom.snp_mat is None:
                par_chrom.snp_mat = pro_chrom.snp_mat.copy()
            else:
                cpp.copy_1D(pro_chrom.snp_mat, par_chrom.snp_mat)
            if par_chrom.snp_pat is None:
                par_chrom.snp_pat = pro_chrom.snp_pat.copy()
            else:
                cpp.copy_1D(pro_chrom.snp_pat, par_chrom.snp_pat)
            founder.mat_id = ls_progeny[x_index].mat_id
            founder.pat_id = ls_progeny[x_index].pat_id
# copy_individual


def nextgen_parents(par_inf, pro_inf, f2_index, m2_index):
    '''
    Copy selected progenies to founder population

    Args:
        par_inf (PopulationInfo): Founder population
        pro_inf (PopulationInfo): Progeny population
        f2_index (list): Arranged index of selected females
        m2_index (list): Arranged index of selected males
    '''
    # female
    copy_individual(par_inf.pop_f, pro_inf.pop_f, f2_index)
    # male
    copy_individual(par_inf.pop_m, pro_inf.pop_m, m2_index)
    # add new parents' ID
    par_inf.new_founder_id()
# nextgen_parent


def start_selection(par_inf, pro_inf, phe_inf, target, method, cross_inf,
                    top_prop, n_family, select_size, max_r):
    '''
    Args:
        par_inf (PopulationInfo): Founder population
        pro_inf (PopulationInfo): Progeny population
        phe_inf (PhenotypeInfo): Phenotype information
        target (str): Selection based on estimated breeding value  ('bv'),
                      phenotype ('phenotype'), or random ('random')
        method (str): How to select from progenies such as mass
                      selection ('mass'), within-family selection
                      ('within-family'), family selection ('family'),
                      or Numertor/Genomic relationship matrix ('FvalueA',
                      'FvalueG')
        cross_inf (numpy.ndarray): Index pairs of female and male parents
        top_prop (float): Select progenies with top X% breeding values
                          in within-family selection. Set 0.0 < top_prop
                          <= 1.0.
        n_family (int): Number of families to be selected
        select_size (tulple): Number of selected founders
        max_r (float): R' among selected individuals is less than this value

    Returns:
        int: If 0, excuted correctly, if 1, terminated irreguraly
    '''
    # Selection target
    select_val = select_value(target, phe_inf)
    # Merge info
    # dict for family, all values are zero
    fam_d = {}
    # list for index, family, select value, 0 (female) or 1 (male)
    summary_pro = []
    # get info for female
    merge_pro_info(pro_inf.pop_f, select_val, fam_d, summary_pro, 0, 0)
    # get info for male
    merge_pro_info(pro_inf.pop_m, select_val, fam_d, summary_pro, pro_inf.n_f,
                   1)
    # mass selection
    if method == 'mass':
        f_index, m_index = mass_selection(summary_pro, select_size)
    # within-family selection
    elif method == 'within-family':
        f_index, m_index = within_family_selection(summary_pro, select_size,
                                                   fam_d, top_prop)
    # family selection
    elif method == 'family':
        f_index, m_index = family_selection(summary_pro, select_size, fam_d,
                                            top_prop, n_family)
    # A/G matrix based selection
    elif method == 'RvalueA':
        f_index, m_index, c_r = fvalue_selection(summary_pro, select_size,
                                                 top_prop, pro_inf.a_mat,
                                                 max_r, pro_inf.n_f)
        if c_r == 1:
            return 1
    elif method == 'RvalueG':
        f_index, m_index, c_r = fvalue_selection(summary_pro, select_size,
                                                 top_prop, pro_inf.g_mat,
                                                 max_r, pro_inf.n_f)
        if c_r == 1:
            return 1
    else:
        sys.exit(f'{method=} is invalid')
    # arrange parents based on cross_info
    # If G matrix is available
    if pro_inf.g_mat is not None:
        k_mat = pro_inf.g_mat
    # only A matrix is available
    else:
        k_mat = pro_inf.a_mat
    f2_index, m2_index = arrange_parents(k_mat, cross_inf, select_size,
                                         f_index, m_index, pro_inf.n_f)
    # Copy selected progenies to founder
    nextgen_parents(par_inf, pro_inf, f2_index, m2_index)
    return 0
# start_selection


def main():
    '''
    main
    '''
    print('A module for selection')
# main


if __name__ == '__main__':
    main()
