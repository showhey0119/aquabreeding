'''
Module for SNPs of parental/progeny populations
'''

import sys
import numpy as np
import msprime as mp
from numba import jit


def identity_by_descent(chrom_inf):
    '''
    Calculate inbreeding coefficiet in a chromosome

    Args:
        chrom_info (ChromInfo class): Chromosome information

    Returns:
        float: Inbreeding coefficient
    '''
    ch_mat = chrom_inf.chrom_mat
    ch_pat = chrom_inf.chrom_pat
    ge_mat = chrom_inf.genotype_mat
    ge_pat = chrom_inf.genotype_pat
    pat_i = mat_i = 0  # index for father and mother
    homo_reg = 0.0  # len IBD regions
    beg_seg = 1  # start position of a segment
    end_seg = -777  # end position of a segment
    while True:
        # check end point
        if ch_mat[mat_i] < ch_pat[pat_i]:
            end_seg = ch_mat[mat_i]
            if ge_mat[mat_i] == ge_pat[pat_i]:
                homo_reg += (end_seg - beg_seg + 1.0)
            mat_i += 1
            beg_seg = end_seg + 1
        else:
            end_seg = ch_pat[pat_i]
            if ge_mat[mat_i] == ge_pat[pat_i]:
                homo_reg += (end_seg - beg_seg + 1)
            pat_i += 1
            beg_seg = end_seg + 1
        if end_seg == chrom_inf.chrom[1]:
            break
    if homo_reg > chrom_inf.chrom[1]:
        print(f'===\n{homo_reg}\n{ch_mat}\n{ge_mat}\n{ch_pat}\n{ge_pat}')
        sys.exit()
    return homo_reg
# proportion_of_ibd


def gene_diversity(snp_array):
    '''
    Calculate gene diversity and the number of segregating sites

    Args:
        snp_array (ndarray): SNP array

    Returns:
        popgen statistics

        - float: Gene divesity
        - int: The number of segregating sites
    '''
    n_row, n_col = np.shape(snp_array)
    hat_h = 0.0
    hat_s = 0
    for i in range(n_col):
        tmp_l = snp_array[:, i]
        f_0 = np.count_nonzero(tmp_l == 0)
        f_1 = np.count_nonzero(tmp_l == 1)
        hat_h += 2.0*f_0*f_1/n_row/(n_row-1.0)
        if f_0 > 0 and f_1 > 0:
            hat_s += 1
    hat_h /= n_col
    return hat_h, hat_s
# gene_diversity


@jit(cache=True)
def progeny_genotype(pro_snp, gen_mat):
    '''
    Convert SNP array into genotype matrix

    Args:
        pro_snp (ndarray): SNP array of progenies
                            (rows: haplotype, columns: loci)
        gen_mat (ndarray): Genotype matrix
                             (rows: genotype, columns: loci)
    '''
    n_rows, n_cols = gen_mat.shape
    for i in range(n_rows):
        for j in range(n_cols):
            gen_mat[i][j] = pro_snp[2*i][j] + pro_snp[2*i+1][j]
# progeny_genotype


def progeny_snp(par_snp, pro_snp, snp_dict, pro_inf):
    '''
    Generating snp array of progeny population

    Args:
        par_snp (ndarray): Founder SNP array
        pro_snp (ndarray): Progeny SNP array
        snp_dict (dict): SNP information
        pro_inf (PopulationInfo class): Progeny population
    '''
    n_snp = par_snp.shape[1]
    # for each SNP
    for j in range(n_snp):
        j_chrom = snp_dict[j]['chrom']
        j_pos = snp_dict[j]['pos']
        # for female
        for i in range(pro_inf.n_popf):
            # get genotype
            gen_mat, gen_pat = pro_inf.pop_f[i].get_genotype(j_chrom, j_pos)
            pro_snp[2*i][j] = par_snp[gen_mat][j]
            pro_snp[2*i+1][j] = par_snp[gen_pat][j]
        # for male
        for i in range(pro_inf.n_popm):
            gen_mat, gen_pat = pro_inf.pop_m[i].get_genotype(j_chrom, j_pos)
            i_2 = i + pro_inf.n_popf
            pro_snp[2*i_2][j] = par_snp[gen_mat][j]
            pro_snp[2*i_2+1][j] = par_snp[gen_pat][j]
# progeny_snp


def generate_snp(n_sample, n_snp):
    '''
    Generate independent SNPs under standard
    Wright-Fisher model

    Args:
        n_smaple (int): Founder size
        n_snp (int): The number of SNPs

    Returns:
        ndarray: SNP matrix
    '''
    rate = (4e-8, 2e-7)  # to simulate almost independent SNPs
    gsa_out = np.empty((0, 2*n_sample), dtype=np.int64)
    gsa_count = 0  # until n_snp
    while True:
        gsa_ts = mp.sim_ancestry(samples=n_sample,
                                 recombination_rate=rate[1],
                                 sequence_length=200,
                                 population_size=10_000)
        gsa_ts = mp.sim_mutations(gsa_ts,
                                  rate=rate[0],
                                  model=mp.BinaryMutationModel(),
                                  discrete_genome=False,
                                  keep=False)
        for g_ts in gsa_ts.variants():
            tmp_gen = g_ts.genotypes
            gsa_out = np.append(gsa_out, np.array([tmp_gen]), axis=0)
            gsa_count += 1
            if gsa_count == n_snp:
                return gsa_out.T
# generate_snp


def main():
    '''
    main
    '''
    print('Module for population genetics')
# main


if __name__ == '__main__':
    main()
