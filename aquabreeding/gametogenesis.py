'''
Module for simulating gametogenesis
'''

import sys
from bisect import bisect_left
import numpy as np


def add_break_point(break_point, chrom_x, geno_x):
    '''
    Add recombination break point at the list of a chromatid

    Args:
        break_point (float): Recombination break point
        chrom_x (list): Chromatid
        geno_x (list): Genotype

    Returns:
        int: index of inserted segment
    '''
    # bisect
    new_index = bisect_left(chrom_x, break_point)
    end_point = int(break_point - 0.5)
    # check if break point already exists
    if new_index > 0 and end_point == chrom_x[new_index-1]:
        return new_index - 1
    # insert
    chrom_x.insert(new_index, end_point)
    new_geno = geno_x[new_index]
    geno_x.insert(new_index, new_geno)
    return new_index
# add_break_point


def exchange_block(ls_p, ls_m, indexp, indexm):
    '''
    Exchange partial list after break point

    Args:
        ls_p (list): Paternal list
        ls_m (list): Maternal list
        indexp (int): Paternal index with break point
        indexm (int): Maternal index with break point
    '''
    ls_p[indexp+1:], ls_m[indexm+1:] = ls_m[indexm+1:], ls_p[indexp+1:]
# exchange_block


def crossing_over(break_point, chrom_p, chrom_m, geno_p, geno_m):
    '''
    Cross-over two chromatids

    * step 1: Get the indices, where crossing-over occurs
    * setp 2: Insert a new segment
    * step 3: Exchange the segments behind the break point

    Args:
        break_point (int): Recombination break point
        chrom_p (list): Paternal chromatid
        chrom_m (list): Maternal chromatid
        geno_p (list): Paternal genotype
        geno_m (list): Maternal genotype
    '''
    new_indexp = add_break_point(break_point, chrom_p, geno_p)
    new_indexm = add_break_point(break_point, chrom_m, geno_m)
    # crossing-over event
    # chromosome
    exchange_block(chrom_p, chrom_m, new_indexp, new_indexm)
    # genotype
    exchange_block(geno_p, geno_m, new_indexp, new_indexm)
# crossing_over


def gameto_genesis(chrom_inf, ex_n_rec):
    '''
    Simulate gametogenesis

    Gametogenesis follows a four-strand model, and allows
    one obligate chiasma.  Crossing-over rate is assumed
    to be constant across chromosomes.

    Args:
        chrom_inf (ChromInfo class): Information of a chromosome
        ex_n_rec (float): Expected number of crossing-over

    Returns:
        Chromosome information

        - list: one of four chromatids,
        - list: genotype

    Note:
        If the expected number of crossing-over event is less than one,
        the script stops.
    '''
    # bivalent chromosomes
    chrom_p1 = chrom_inf.chrom_pat.copy()
    chrom_p2 = chrom_inf.chrom_pat.copy()
    chrom_m1 = chrom_inf.chrom_mat.copy()
    chrom_m2 = chrom_inf.chrom_mat.copy()
    geno_p1 = chrom_inf.genotype_pat.copy()
    geno_p2 = chrom_inf.genotype_pat.copy()
    geno_m1 = chrom_inf.genotype_mat.copy()
    geno_m2 = chrom_inf.genotype_mat.copy()
    # crossing-over
    l_chrom = chrom_inf.chrom[1]
    # '1+' means obligate chiasma
    if ex_n_rec > 1.0:
        n_rec = 1 + np.random.poisson(ex_n_rec-1)
    # if crossing-over rate = 0
    elif ex_n_rec == 0:
        n_rec = 0
    # if chromosome is very short
    else:
        n_rec = 1
        sys.exit('Expected num crossing over events is less than 1')
    # if 100.5, crossing-over occurs between 100th and 101th sites
    break_point = np.random.randint(low=1, high=l_chrom, size=n_rec) + 0.5
    bivalent = np.random.randint(4, size=n_rec+1)
    # crossing-over
    for i in range(n_rec):
        # which sister-chromatids are recombined
        if bivalent[i] == 0:
            crossing_over(break_point[i], chrom_p1, chrom_m1, geno_p1, geno_m1)
        elif bivalent[i] == 1:
            crossing_over(break_point[i], chrom_p1, chrom_m2, geno_p1, geno_m2)
        elif bivalent[i] == 2:
            crossing_over(break_point[i], chrom_p2, chrom_m1, geno_p2, geno_m1)
        elif bivalent[i] == 3:
            crossing_over(break_point[i], chrom_p2, chrom_m2, geno_p2, geno_m2)
    # which chromaid is inherited
    if bivalent[n_rec] == 0:
        return chrom_p1, geno_p1
    if bivalent[n_rec] == 1:
        return chrom_p2, geno_p2
    if bivalent[n_rec] == 2:
        return chrom_m1, geno_m1
    if bivalent[n_rec] == 3:
        return chrom_m2, geno_m2
    return -777
# gameto_genesis


def main():
    '''
    main
    '''
    print('gametogenesis.py')
    print('Module for gametogenesis')
# main


if __name__ == '__main__':
    main()
