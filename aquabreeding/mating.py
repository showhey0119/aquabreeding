'''
A module for mating
'''

import re
import sys
import numpy as np
from aquabreeding import gametogenesis as gg


def set_mating_design_partial(select_size, design):
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
        select_size (tuple): Nos. female and male parents
        design (str): 1x(int)

    Returns:
        numpy.ndarray: Index pairs of female and male parents
    '''
    n_cross = re.findall(r'^1x(\d+)$', design)
    if len(n_cross) != 1:
        sys.exit('cross_design should be \'1x[int]\'')
    n_cross = int(n_cross[0])
    if n_cross == select_size[0]:
        sys.exit('Use \'full\' for design instead')
    if n_cross > select_size[0]:
        sys.exit('No. cross is too large')
    pair_i = np.empty((0, 2), dtype=np.int32)
    for i in range(select_size[0]):
        for j in range(i, i+n_cross):
            if j < select_size[0]:
                j_2 = j
            else:
                j_2 = j - select_size[0]
            pair_i = np.append(pair_i, np.array([[i, j_2]]), axis=0)
    return pair_i.astype(np.int32)
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
        select_size (tuple): Nos. female and male parents

    Returns:
        numpy.ndarray: index pairs of female and male parents
    '''
    pair_i = np.empty((0, 2), dtype=np.int32)
    for i in range(select_size[0]):
        for j in range(select_size[1]):
            pair_i = np.append(pair_i, np.array([[i, j]]), axis=0)
    return pair_i.astype(np.int32)
# set_mating_design_full


def set_mating_design(design, select_size):
    '''
    Generate cross design information

    Args:
        design (unknown): Implemented mating design (str) or index pairs of
                          female and male parents (numpy.ndarray)
        select_size (tuple): Nos. selected female and male parents

    Returns:
        numpy.ndarray: Female index, male index
    '''
    # implemented
    if isinstance(design, str):
        if re.compile(r'^1x\d+$').match(design):
            return set_mating_design_partial(select_size, design)
        if design == 'full':
            return set_mating_design_full(select_size)
    # custom
    if isinstance(design, np.ndarray):
        c_design = np.shape(design)[1]
        if c_design != 2:
            sys.exit('No. columns of mating design should be 2')
        return design.astype(np.int32)
    sys.exit('design should be 1x(int), full, or numpy.ndarray')
# set_mating_design


def mating_process(cross_inf, par_inf, pro_ls, n_pro):
    '''
    Produce female/male progenies

    Args:
        cross_inf (numpy.ndarray): Index pairs of female and male parents
        par_inf (PopulationInfo): Founder population
        pro_ls (list): List of IndividualInfo class of female/male
        n_pro (int): No. progenies
    '''
    n_cross = np.shape(cross_inf)[0]
    # No. progenies in each cross
    family_size = np.random.multinomial(n_pro, [1.0/n_cross]*n_cross)
    # index for progeny
    i_pro = 0
    # for parental female index, male index, no. progeny
    for i, j, p_size in zip(cross_inf[:, 0], cross_inf[:, 1], family_size):
        for _ in range(p_size):
            gg.produce_progeny(par_inf.pop_f[i], par_inf.pop_m[j],
                               pro_ls[i_pro])
            i_pro += 1
    if i_pro != n_pro:
        sys.exit('something wrong in produce_progeny_female_male')
# maiting_process


def start_mating(cross_inf, par_inf, pro_inf):
    '''
    Produce progenies by crossing founder/parental individuals

    Args:
        cross_inf (numpy.ndarray): Index pairs of female and male parents
        par_inf (PopulationInfo): Founder population
        pro_inf (PopulationInfo): Progeny population
    '''
    # produce female progeny
    mating_process(cross_inf, par_inf, pro_inf.pop_f, pro_inf.n_f)
    # produce male progeny
    mating_process(cross_inf, par_inf, pro_inf.pop_m, pro_inf.n_m)
# start_mating


def main():
    '''
    main
    '''
    print('A module for mating')
# main


if __name__ == '__main__':
    main()
