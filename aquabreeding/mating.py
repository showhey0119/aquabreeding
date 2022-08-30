'''
module for mating
'''

import sys
import random
import operator
import numpy as np


def error_select_by(mess):
    '''
    If select_by is unexpected, error message is shown
    in mass_selection or within_family_selection
    args: mess: 'mass' or 'within family' (str)
    '''
    print(f'Error in {mess} selection')
    print('select_by should be')
    print(' \'bv\': selection based on breeding values')
    print(' \'phenotype\': selection based on phenotypes')
    print(' \'no\': no selection')
    sys.exit()
    # error_select_by


def select_wf(n_p, n_f, p_pop, top_prop, target=None):
    '''
    select progenies as the next founder population
    retain pedigrees by within family selection
    args: n_p: progeny size (int)
          n_f: founder size (int)
          p_pop: progeny population (list of class)
          target: list of bv or phenotype (str)
                  if None, random number array is generated
          top_prop: top {100*top_prop} percentage of progenies are used.
                    0.0 ~ 1.0 (float)
    return: IDs of selected progenies (list)
    '''
    # if no selection
    if target is None:
        target = np.random.rand(n_p)
    # check top_prop
    if top_prop > 1.0 or top_prop < 0.0:
        sys.exit('top_prop should be 0.0 ~ 1.0')
    if top_prop * n_p < n_f:
        sys.exit('top_prop * progeny_size should be larger than founder_size')
    # check family
    fam_d = {}
    tmp_bv = []
    for i in range(n_p):
        if p_pop[i].pat_id < p_pop[i].mat_id:
            fam_n = f'f{p_pop[i].pat_id}x{p_pop[i].mat_id}'
        else:
            fam_n = f'f{p_pop[i].mat_id}x{p_pop[i].pat_id}'
        tmp_bv.append([i, fam_n, target[i]])
        fam_d.setdefault(fam_n, 0)
    tmp_bv.sort(reverse=True, key=operator.itemgetter(2))
    # print(fam_d)
    # for i in tmp_bv:
    #    print(i)
    # get ID of larget bv within families
    i_res = []
    n_res = 0
    n_tolerant = 0  # how many times progenies are selected in a family
    n_p2 = int(top_prop*n_p)
    if n_p2 < n_f:
        print(f'{n_p2}, {n_f}')
        sys.exit('Error in select_wf')
    while True:
        for i in range(n_p2):
            if fam_d[tmp_bv[i][1]] == n_tolerant and tmp_bv[i][0] != -7:
                i_res.append(tmp_bv[i][0])
                n_res += 1
                fam_d[tmp_bv[i][1]] += 1
                tmp_bv[i][0] = -7  # avoid picking up the same progeny
                # print(tmp_bv[i])
                if n_res == n_f:
                    random.shuffle(i_res)
                    return i_res
        n_tolerant += 1
    # select_wf


def selected_progeny(n_p, n_f, target=None):
    '''
    get index of progeny with largest bv or phenotype
    args: n_p: progeny size
          n_f: founder size
          target: selection based on bv, phenotype, or None
    '''
    if target is None:
        target = np.random.rand(n_p)
    id_bv = [[i, target[i]] for i in range(n_p)]
    id_bv.sort(reverse=True, key=operator.itemgetter(1))
    i_res = []
    n_res = 0
    for i in range(n_p):
        i_res.append(id_bv[i][0])
        n_res += 1
        if n_res == n_f:
            random.shuffle(i_res)
            return i_res
    return None
    # selected_progeny


def factorial_design(cross_design, n_founder, n_progeny):
    '''
    mating design in factorial cross
    args: cross_design (str), 1x1, 1x2, 2x2, full
          n_founder (int)
          n_progeny (int)
    return list([pat ID, mat ID, replicate], [],...,[]])
    '''
    if cross_design == '1x1':
        #   1 3 5 7
        # 0 x
        # 2   x
        # 4     x
        # 6       x
        # the number of mating pairs
        n_cross = int(n_founder/2)
        # indices of mating pairs
        pair_i = []
        for i in range(0, n_founder, 2):
            pair_i.append([i, i+1, 0])
        # 1x1
    elif cross_design == '1x2':
        #   1 3 5 7
        # 0 x x
        # 2   x x
        # 4     x x
        # 6 x     x
        n_cross = n_founder
        pair_i = []
        for i in range(0, n_founder, 2):
            for j in range(i+1, i+4, 2):
                if j == n_founder+1:
                    j = 1
                pair_i.append([i, j, 0])
        # 1x2
    elif cross_design == '2x2':
        #   1 3 5 7 9
        # 0 x x
        # 2 x x
        # 4     x x
        # 6     x x
        # 8         x
        if n_founder % 4 == 0:
            n_cross = n_founder
        else:
            n_cross = n_founder - 1
        pair_i = []
        for k in range(0, n_cross, 4):
            if k == n_cross-1:
                pair_i.append([k, k+1, 0])
                break
            for i in range(k, k+3, 2):
                for j in range(k+1, k+4, 2):
                    pair_i.append([i, j, 0])
        # 2x2
    elif cross_design == 'full':
        #   1 3 5 7
        # 0 x x x x
        # 2 x x x x
        # 4 x x x x
        # 6 x x x x
        n_cross = n_founder ** 2 / 4
        pair_i = []
        for i in range(0, n_founder, 2):
            for j in range(1, n_founder, 2):
                pair_i.append([i, j, 0])
        # full
    else:
        sys.exit('cross design should be 1x1, 1x2, 2x2, or full')
    # the expected numbers of progenies in cross pairs
    # are the same
    freq_std = 1.0
    ex_p = 1.0/n_cross
    n_total = n_progeny
    for pair_j in pair_i:
        binom_p = ex_p/freq_std
        binom_p = min(binom_p, 1.0)
        pair_j[2] = np.random.binomial(n_total, binom_p)
        freq_std -= ex_p
        n_total -= pair_j[2]
    return pair_i
    # factorial_design


def different_randint_set(r_max, n_set):
    '''
    obtain {n_set} sets of different randint (0 ~ {r_max-1})
    for random_mating
    return  np.ndarray as
    [[s0_1, s0_2], [s1_1, s1_2], [s2_1, s2_2]...]]
    '''
    res_r = np.empty((0, 2), dtype=np.int)
    while True:
        tmp_r = np.random.randint(r_max, size=(2*n_set, 2))
        tmp_r = tmp_r[tmp_r[:, 0] != tmp_r[:, 1]]
        res_r = np.append(res_r, tmp_r, axis=0)
        if res_r.shape[0] >= n_set:
            return res_r[0:n_set, :]
    # different_randint_set


def main():
    '''
    main
    '''
    print('The module for mating design')
    # main


if __name__ == '__main__':
    main()
