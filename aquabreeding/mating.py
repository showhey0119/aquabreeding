'''
Module for mating
'''

import re
import sys
import numpy as np
from numba import jit


def set_mating_design_partial(select_size, cross_design):
    '''
    Cross design information for partial factorial cross

    ::

                    male
                  1x1            1x2            1x3
                  0 1 2 3 4      0 1 2 3 4      0 1 2 3 4
                0 x            0 x x          0 x x x
                1   x          1   x x        1   x x x
        female  2     x        2     x x      2     x x x   ...
                3       x      3       x x    3 x     x x
                4         x    4 x       x    4 x x     x

    Args:
        par_inf (PopulationInfo class): Founder population

    Returns:
        ndarray: female index, male index, 0, 0
    '''
    n_cross = re.findall(r'^1x(\d+)$', cross_design)
    if len(n_cross) != 1:
        sys.exit('cross_design should be \'1x[int]\'')
    n_cross = np.int(n_cross[0])
    if n_cross == select_size[0]:
        sys.exit('Use \'full\' for cross_design instead')
    if n_cross > select_size[0]:
        sys.exit('No. cross is too large')
    pair_i = np.empty((0, 4), dtype=np.int64)
    for i in range(select_size[0]):
        for j in range(i, i+n_cross):
            if j < select_size[0]:
                j_2 = j
            else:
                j_2 = j - select_size[0]
            pair_i = np.append(pair_i, np.array([[i, j_2, 0, 0]]), axis=0)
    return pair_i
# set_mating_design_partial


def set_mating_design_full(select_size):
    '''
    Cross design information for full factorial cross

    ::

                   male
                 0 1 2 3 4
               0 x x x x x
               1 x x x x x
        female 2 x x x x x
               3 x x x x x
               4 x x x x x

    Args:
        par_inf (PopulationInfo class): Founder population

    Returns:
        ndarray: female index, male index, 0, 0
    '''
    pair_i = np.empty((0, 4), dtype=np.int64)
    for i in range(select_size[0]):
        for j in range(select_size[1]):
            pair_i = np.append(pair_i, np.array([[i, j, 0, 0]]), axis=0)
    return pair_i
# set_mating_design_full


def set_mating_design(select_size, cross_design):
    '''
    Generate cross design information

    Args:
        par_inf (PopulationInfo class): Founder population
        cross_design (str): '1x(int)' partial or 'full' factorial mating

    Returns:
        ndarray: Female index, male index, 0, 0
    '''
    if re.compile(r'^1x\d+$').match(cross_design):
        pair_i = set_mating_design_partial(select_size, cross_design)
    elif cross_design == 'full':
        pair_i = set_mating_design_full(select_size)
    else:
        sys.exit('cross_design should be 1x(int), or full')
    return pair_i
# mate_design


@jit(cache=True)
def random_allocation(cross_inf, n_total, tag):
    '''
    Randomly allocate the numbers of progenies in each cross

    Args:
        cross_inf (ndarray): Female index, male index,
                             no. female progeny, no. male progeny
        n_total (int): No. male/female progenies
        tag (int): Where female/male progeny num is stored in cross_inf.
                   2 for female, 3 for male.
    '''
    n_cross = cross_inf.shape[0]
    freq_std = 1.0  # normalized total freq
    ex_p = 1.0/n_cross  # expected frequency
    for i in range(n_cross):
        binom_p = ex_p/freq_std  # binomial p is normalized
        binom_p = min(binom_p, 1.0)  # avoid 1.00000000001
        cross_inf[i][tag] = np.random.binomial(n_total, binom_p)
        freq_std -= ex_p
        n_total -= cross_inf[i][tag]
# random_allocation


def produce_progeny_female_male(cross_inf, par_inf, pro_ls,
                                n_pro, tag):
    '''
    Produce female/male progenies

    Args:
        cross_inf (ndarray): Female index, male index,
                             no. female progeny, no. male progeny
        par_inf (PopulationInfo class): Founder population
        pro_ls (list): List of IndividualInfo class of female/male
        n_pro (int): Length of pro_ls
        tag (int): Where female/male progeny num is stored in cross_inf.
                   2 for female, 3 for male
    '''
    # index for progeny
    i_pro = 0
    # for parental female index, male index, no. progeny
    for i, j, p_size in zip(cross_inf[:, 0], cross_inf[:, 1],
                            cross_inf[:, tag]):
        for _ in range(p_size):
            # for each chromosome
            for c_id in range(par_inf.chrom[0]):
                # female
                gamete_1, geno_1 = par_inf.pop_f[i].gamete(c_id, 0)
                # male
                gamete_2, geno_2 = par_inf.pop_m[j].gamete(c_id, 1)
                # copy gamete
                pro_ls[i_pro].copy_gametes(gamete_1, gamete_2, geno_1, geno_2,
                                           c_id)
            # parents' ID
            pro_ls[i_pro].mat_id = par_inf.pop_f[i].ind_id
            pro_ls[i_pro].pat_id = par_inf.pop_m[j].ind_id
            i_pro += 1
    if i_pro != n_pro:
        sys.exit('something wrong in produce_progeny_female_male')
# produce_progeny_female_male


def produce_progeny(cross_inf, par_inf, pro_inf, r_a):
    '''
    Produce progenies by crossing founder/parental individuals

    Args:
        cross_inf (ndarray): Female index, male index,
                             no. female progeny, no. male progeny
        par_inf (PopulationInfo class): Founder population
        pro_inf (PopulationInfo class): Progeny population
        ra (bool): random allocation or not

    Note:
        The numbers of female/male progenies are stored in
        2nd and 3rd colums in cross_inf, respectively.
    '''
    if r_a:
        # simulate number of progenies in each cross for female
        random_allocation(cross_inf, pro_inf.n_popf, 2)
        # for male
        random_allocation(cross_inf, pro_inf.n_popm, 3)
    # produce female progeny
    produce_progeny_female_male(cross_inf, par_inf, pro_inf.pop_f,
                                pro_inf.n_popf, 2)
    # produce male progeny
    produce_progeny_female_male(cross_inf, par_inf, pro_inf.pop_m,
                                pro_inf.n_popm, 3)
# produce_progeny


def main():
    '''
    main
    '''
    print('mating.py')
    print('Module for mating')
# main


if __name__ == '__main__':
    main()
