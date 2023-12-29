'''
A module for simulating gametogenesis

Notes:
    SNP position starts with zero
'''

import sys
import numpy as np
from aquabreeding import aquacpp as cpp


def new_gamete(chrom_info, ex_n_rec):
    '''
    Simulate gametogenesis

    Gametogenesis follows a four-strand model, and allows
    one obligate chiasma.  Crossing-over rate is assumed
    to be constant across chromosomes.

    Args:
        chrom_info (ChromInfo class): Information of a chromosome
        ex_n_rec (float): Expected number of crossing-over

    Returns:
        numpy.ndarray: SNPs of a new gamete

    Note:
        If the expected number of crossing-over event is less than one,
        the script stops.
    '''
    # bivalent chromosomes
    chrom_p1 = chrom_info.snp_pat.copy()
    chrom_p2 = chrom_info.snp_pat.copy()
    chrom_m1 = chrom_info.snp_mat.copy()
    chrom_m2 = chrom_info.snp_mat.copy()
    # crossing-over
    l_chrom = chrom_info.chrom[1]
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
    break_point = np.random.randint(low=0, high=l_chrom - 1, size=n_rec) + 0.5
    bivalent = np.random.randint(4, size=n_rec+1)
    # crossing-over
    for i in range(n_rec):
        # which sister-chromatids are recombined
        if bivalent[i] == 0:
            cpp.crossing_over(break_point[i], chrom_p1, chrom_m1,
                              chrom_info.position)
        elif bivalent[i] == 1:
            cpp.crossing_over(break_point[i], chrom_p1, chrom_m2,
                              chrom_info.position)
        elif bivalent[i] == 2:
            cpp.crossing_over(break_point[i], chrom_p2, chrom_m1,
                              chrom_info.position)
        elif bivalent[i] == 3:
            cpp.crossing_over(break_point[i], chrom_p2, chrom_m2,
                              chrom_info.position)
    # which chromaid is inherited
    if bivalent[n_rec] == 0:
        return chrom_p1
    if bivalent[n_rec] == 1:
        return chrom_p2
    if bivalent[n_rec] == 2:
        return chrom_m1
    if bivalent[n_rec] == 3:
        return chrom_m2
    return -777
# new_gameto


def produce_progeny(par_f, par_m, pro_x):
    '''
    Produce progeny

    Args:
        par_f (IndividualInfo): Female parent
        par_m (IndividualInfo): Mmale parent
        pro_x (IndividualInfo): Progeny
    '''
    # for each chromosome
    for i in range(par_f.chrom[0]):
        # copy position
        if pro_x.chrom_ls[i].position is None:
            pro_x.chrom_ls[i].position = par_f.chrom_ls[i].position.copy()
        # gamete from female
        new_gamete1 = new_gamete(par_f.chrom_ls[i], par_f.n_recf)
        if pro_x.chrom_ls[i].snp_mat is None:
            pro_x.chrom_ls[i].snp_mat = new_gamete1.copy()
        else:
            cpp.copy_1D(new_gamete1, pro_x.chrom_ls[i].snp_mat)
        # gamete from male
        new_gamete2 = new_gamete(par_m.chrom_ls[i], par_m.n_recm)
        if pro_x.chrom_ls[i].snp_pat is None:
            pro_x.chrom_ls[i].snp_pat = new_gamete2.copy()
        else:
            cpp.copy_1D(new_gamete2, pro_x.chrom_ls[i].snp_pat)
    # Set ID
    pro_x.mat_id = par_f.sample_id
    pro_x.pat_id = par_m.sample_id
    pro_x.sample_id = -1
# produce_progeny


def main():
    '''
    main
    '''
    print('A module for gametogenesis')
# main


if __name__ == '__main__':
    main()
