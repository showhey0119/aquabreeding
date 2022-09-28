'''
module for simulating gametogenesis
'''

import numpy as np


def add_break_point(ab_posi, ab_ls):
    '''
    add recombination break point at the list of a chromatid

    Args:
        ab_posi (int): recombination break point
        ab_ls (list): a chromatid

    Returns:
        ab_i (int): the index of inserted segment
    '''
    for ab_i, ab_tmp in enumerate(ab_ls):
        # add new break point
        if ab_tmp[0] < ab_posi < ab_tmp[1]:
            ab_ls.insert(ab_i, [ab_tmp[0], int(ab_posi-0.5), ab_tmp[2]])
            ab_ls[ab_i+1][0] = int(ab_posi+0.5)
            return ab_i
        # when break point already exists
        if int(ab_posi-0.5) == ab_tmp[1]:
            return ab_i
    return -777
# add_break_point


def crossing_over(es_posi, es_ls1, es_ls2):
    '''
    cross over two chromatids

    step 1: get the indices, where crossing-over occurs
    step 2: exchange the segments behind the break point

    Args:
        es_posi (int): recombination break point
        es_ls1 (list): a chromatid
        es_ls2 (list): another chromatid
    '''
    es_index1 = add_break_point(es_posi, es_ls1)
    es_index2 = add_break_point(es_posi, es_ls2)
    # crossing-over event
    es_ls1_behind = es_ls1[es_index1+1:]
    del es_ls1[es_index1+1:]
    es_ls2_behind = es_ls2[es_index2+1:]
    del es_ls2[es_index2+1:]
    es_ls1.extend(es_ls2_behind)
    es_ls2.extend(es_ls1_behind)
# exchange_segment


def gameto_genesis(ch_pat, ch_mat, ex_n_rec, l_chrom):
    '''
    Simulate gametogenesis

    Gametogenesis follows four-strand model, and allows
    one obligate chiasma.  Crossing-over rate is assumed
    to be constant across the genome

    Args:
        ch_pat (list): paternal chromosome
        ch_mat (list): maternal chromosome
        ex_n_rec (float): expected number of crossing-over
        l_chrom (int): chrom len (bp)

    Returns:
        pat/mat_1/2 (list): one of four chromatids
    '''
    # bivalent chromosomes
    pat_1 = [[i_1[0], i_1[1], i_1[2]] for i_1 in ch_pat]
    pat_2 = [[i_1[0], i_1[1], i_1[2]] for i_1 in ch_pat]
    mat_1 = [[i_1[0], i_1[1], i_1[2]] for i_1 in ch_mat]
    mat_2 = [[i_1[0], i_1[1], i_1[2]] for i_1 in ch_mat]
    # crossing-over
    # '1+' means obligate chiasma
    if ex_n_rec > 1.0:
        n_rec = 1 + np.random.poisson(ex_n_rec-1)
    else:
        n_rec = 1
    break_point = np.random.randint(low=1, high=l_chrom, size=n_rec)+0.5
    bivalent = np.random.randint(4, size=n_rec+1)
    # crossing-over
    for i in range(n_rec):
        # which sister-chromatids are recombined
        if bivalent[i] == 0:
            crossing_over(break_point[i], pat_1, mat_1)
        elif bivalent[i] == 1:
            crossing_over(break_point[i], pat_1, mat_2)
        elif bivalent[i] == 2:
            crossing_over(break_point[i], pat_2, mat_1)
        elif bivalent[i] == 3:
            crossing_over(break_point[i], pat_2, mat_2)
    # which chtomaid is inherited
    if bivalent[n_rec] == 0:
        return pat_1
    if bivalent[n_rec] == 1:
        return pat_2
    if bivalent[n_rec] == 2:
        return mat_1
    if bivalent[n_rec] == 3:
        return mat_2
    return -777
# gameto_genesis


def main():
    '''
    main
    '''
    print('Module for gametogenesis')
# main


if __name__ == '__main__':
    main()
