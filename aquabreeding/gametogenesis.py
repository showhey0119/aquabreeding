'''
module for simulating gametogenesis
'''

# import sys
import numpy as np
from numba import jit


@jit(cache=True)
def where_break_point(w_posi, w_ls):
    '''
    search the segment where crossing-over occurs
    '''
    for i in range(w_ls.shape[0]):
        if w_ls[i][0] < w_posi < w_ls[i][1]:
            return i
        if int(w_posi-0.5) == w_ls[i][1]:
            return -7
    return None
    # where_break_point


@jit(cache=True)
def where_break_point2(w_posi, w_ls):
    '''
    search the segment where crossing-over occurs
    call this function when crossing-over occurs
    at the alrealy-existing break point
    '''
    for i in range(w_ls.shape[0]):
        if int(w_posi-0.5) == w_ls[i][1]:
            return i
    return None
    # where_break_point2


@jit(cache=True)
def add_break_point(ab_i, ab_ls, ab_posi):
    '''
    add recombination break point
    '''
    for i in range(ab_ls.shape[0]-1, -1, -1):
        if ab_ls[i][0] == -7:
            continue
        if i == ab_i:
            ab_ls[i+1][0] = int(ab_posi+0.5)
            ab_ls[i+1][1] = ab_ls[i][1]
            ab_ls[i+1][2] = ab_ls[i][2]
            ab_ls[i][1] = int(ab_posi-0.5)
            break
        ab_ls[i+1][0] = ab_ls[i][0]
        ab_ls[i+1][1] = ab_ls[i][1]
        ab_ls[i+1][2] = ab_ls[i][2]
    # add_break_point


@jit('void(i8, i8, i8[:,:], i8[:,:])', cache=True)
def exchange_segments(e_i1, e_i2, e_ls1, e_ls2):
    '''
    exchange the genomic segments behind the
    recombination break point
    '''
    l_ls1 = e_ls1.shape[0]
    l_ls2 = e_ls2.shape[0]
    # ls1-ls2
    e_ls1_behind = e_ls1[e_i1+1:].copy()
    e_j = e_i2+1
    for i in range(e_i1+1, l_ls1):
        if e_j >= l_ls2:
            e_ls1[i][0] = -7
            e_ls1[i][1] = -7
            e_ls1[i][2] = -7
        else:
            e_ls1[i][0] = e_ls2[e_j][0]
            e_ls1[i][1] = e_ls2[e_j][1]
            e_ls1[i][2] = e_ls2[e_j][2]
        e_j += 1
    # ls2-ls1
    e_j = 0
    l_ls1_b = e_ls1_behind.shape[0]
    for i in range(e_i2+1, l_ls2):
        if e_j >= l_ls1_b:
            e_ls2[i][0] = -7
            e_ls2[i][1] = -7
            e_ls2[i][2] = -7
        else:
            e_ls2[i][0] = e_ls1_behind[e_j][0]
            e_ls2[i][1] = e_ls1_behind[e_j][1]
            e_ls2[i][2] = e_ls1_behind[e_j][2]
        e_j += 1
    # exchange_segments


def crossing_over(es_posi, es_ls1, es_ls2):
    '''
    args: recombination break point (int)
          chromatid 1 (np.ndarray)
          chromatid 2 (np.ndarray)
    step 1: get the indices, where crossing-over occurs
    setp 2: add recombination break point
    step 3: exchange the segments behind the break point
    '''
    # search segment where crossing-over occurs
    # first chromatid
    es_i1 = where_break_point(es_posi, es_ls1)
    if es_i1 == -7:
        # if the break point already exists
        es_i1 = where_break_point2(es_posi, es_ls1)
    else:
        add_break_point(es_i1, es_ls1, es_posi)
    # second chromatid
    es_i2 = where_break_point(es_posi, es_ls2)
    if es_i2 == -7:
        es_i2 = where_break_point2(es_posi, es_ls2)
    else:
        add_break_point(es_i2, es_ls2, es_posi)
    # crossing-over
    exchange_segments(es_i1, es_i2, es_ls1, es_ls2)
    # crossing_over


@jit(cache=True)
def first_7(ls_7):
    '''
    search not-used memory allocation
    '''
    for i in range(ls_7.shape[0]):
        if ls_7[i][0] == -7:
            return i
    return None
    # first_7


def gameto_genesis(ch_pat, ch_mat, ex_n_rec, l_chrom):
    '''
    a function for simulating gametogenesis
    args: paternal chromosome (np.array), maternal chromosome (np.array),
          expected number of crossing-over (float)
          chrom len (bp) (int)
    return one of four chromatids (np.ndarray)
    '''
    # bivalent chromosomes
    pat_1 = ch_pat.copy()
    pat_2 = ch_pat.copy()
    mat_1 = ch_mat.copy()
    mat_2 = ch_mat.copy()
    # crossing-over
    # '1+' means obligate chiasma
    n_rec = 1 + np.random.poisson(ex_n_rec-1)
    break_point = np.random.randint(low=1, high=l_chrom, size=n_rec)+0.5
    bivalent = np.random.randint(4, size=n_rec+1)
    # memory allocation
    l_max = np.max((ch_pat.shape[0], ch_mat.shape[0]))
    tmp_7 = [[-7, -7, -7] for i in range(n_rec+l_max)]
    pat_1 = np.concatenate([pat_1, tmp_7], axis=0)
    pat_2 = np.concatenate([pat_2, tmp_7], axis=0)
    mat_1 = np.concatenate([mat_1, tmp_7], axis=0)
    mat_2 = np.concatenate([mat_2, tmp_7], axis=0)
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
        g_i = first_7(pat_1)
        return pat_1[:g_i]
    if bivalent[n_rec] == 1:
        g_i = first_7(pat_2)
        return pat_2[:g_i]
    if bivalent[n_rec] == 2:
        g_i = first_7(mat_1)
        return mat_1[:g_i]
    if bivalent[n_rec] == 3:
        g_i = first_7(mat_2)
        return mat_2[:g_i]
    return -777
    # gametogenesis


def main():
    '''
    main
    '''
    print('The module for gametogenesis')
    # main


if __name__ == '__main__':
    main()
