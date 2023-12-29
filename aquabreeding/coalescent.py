'''
A module for coalescent simulation

Notes:
    * SNP position starts with zero
'''

import sys
import numpy as np
import msprime as mp


def copy_snp(chrom_inf, tmp_pos, tmp_snp, i_c):
    '''
    Copy SNPs

    Args:
        chrom_inf (ChromInfo): Chromosome info
        tmp_pos (numpy.ndarray): Position of SNPs
        tmp_snp (numpy.ndarray): SNP matrix
        i_c (int): Counter

    Returns:
        int: Counter
    '''
    chrom_inf.position = tmp_pos.copy()
    chrom_inf.snp_mat = tmp_snp[i_c].copy()
    i_c += 1
    chrom_inf.snp_pat = tmp_snp[i_c].copy()
    i_c += 1
    return i_c
# copy_snp


def add_founder_snp(par_inf, snp_mat):
    '''
    Add SNPs to the founder

    Args:
        par_inf (PopulationInfo): founder population
        snp_mat (numpy.ndarray): SNP matrix (np. haplotypes x np. SNP)
    '''
    l_chrom = par_inf.chrom[1]
    # No. SNPs
    c_snp = np.shape(snp_mat)[1]
    # Divide SNPs into chromosomes
    n_chrom = par_inf.chrom[0]
    n_divide = np.random.multinomial(c_snp, [1.0/n_chrom]*n_chrom)
    s_st = None
    s_en = None
    for i in range(n_chrom):
        # Divide SNPs
        if s_st is None:
            s_st = 0
        else:
            s_st = s_en
        s_en = s_st + n_divide[i]
        tmp_snp = snp_mat[:, s_st:s_en]
        # position
        tmp_pos = np.random.randint(low=0, high=l_chrom, size=n_divide[i],
                                    dtype=np.int32)
        tmp_pos.sort()
        # Store the data
        i_c = 0
        for individual in par_inf.pop_f:
            i_c = copy_snp(individual.chrom_ls[i], tmp_pos, tmp_snp, i_c)
        for individual in par_inf.pop_m:
            i_c = copy_snp(individual.chrom_ls[i], tmp_pos, tmp_snp, i_c)
# add_founder_snp


def get_index_structured(n_pop, n_female, n_male):
    '''
    Get index of females and males in the result of
    msprime with population structure

    Args:
        n_pop (int): The number of population
        n_female (tuple): The numbers of females in each population
        n_male (tuple): The numbers of males in each population

    Returns:
        list, list: lists of female/male index
    '''
    tmp_i = 0
    f_id = []
    m_id = []
    for i in range(n_pop):
        # female
        for _ in range(2*n_female[i]):
            f_id.append(tmp_i)
            tmp_i += 1
        # male
        for _ in range(2*n_male[i]):
            m_id.append(tmp_i)
            tmp_i += 1
    return f_id, m_id
# get_index_structured


def structured_population(n_snp, gblup, n_pop, fst_value, n_female, n_male):
    '''
    Run msprime with structured populations

    Args:
        n_snp (int): No. causal SNPs
        gblup (int): No. neutral SNPs
        n_pop (int): The number of populations
        fst_value (float): Average Fst value among populations
        n_female (tuple): The numbers of females in each population
        n_male (tuple): The number of males in each population

    Returns:
        numpy.ndarray: SNP matrix (no. haplotypes x no. SNPs)
    '''
    # Total number of SNPs
    if gblup is None:
        n_total = n_snp
    else:
        n_total = n_snp + gblup
    rate = (4e-8, 2e-7)  # to simulate almost independent SNPs
    # index of females and males haplotype
    sum_f = 2*sum(n_female)
    sum_m = 2*sum(n_male)
    f_id, m_id = get_index_structured(n_pop, n_female, n_male)
    # Set demography
    each_pop = 10_000
    time_div = 2.0 * each_pop * (1.0 / (1.0 - fst_value) - 1.0)
    pop_size = [each_pop]*(n_pop+1)
    demography = mp.Demography.isolated_model(initial_size=pop_size)
    demography.add_population_split(time=time_div,
                                    derived=list(range(n_pop)),
                                    ancestral=n_pop)
    # sample size
    sample_dict = {}
    for i in range(n_pop):
        sample_dict[i] = n_female[i] + n_male[i]
    # output array
    out_f = np.empty((0, sum_f), dtype=np.int32)
    out_m = np.empty((0, sum_m), dtype=np.int32)
    # run msprime
    gsa_count = 0
    while True:
        gsa_ts = mp.sim_ancestry(samples=sample_dict,
                                 recombination_rate=rate[1],
                                 sequence_length=200,
                                 demography=demography)
        gsa_ts = mp.sim_mutations(gsa_ts,
                                  rate=rate[0],
                                  model=mp.BinaryMutationModel(),
                                  discrete_genome=False,
                                  keep=False)
        for g_ts in gsa_ts.variants():
            tmp_gen = g_ts.genotypes
            # skip monomorphic
            a_fq = np.sum(tmp_gen)
            if a_fq in (0, sum_f + sum_m):
                continue
            tmp_f = tmp_gen[f_id]
            tmp_m = tmp_gen[m_id]
            out_f = np.append(out_f, np.array([tmp_f]), axis=0)
            out_m = np.append(out_m, np.array([tmp_m]), axis=0)
            gsa_count += 1
            if gsa_count == n_total:
                gsa_out = np.concatenate([out_f, out_m], axis=1)
                return gsa_out.T.copy()
# structured_population


def run_msprime(n_snp, gblup, n_sample):
    '''
    Run msprime under a given model

    Args:
        par_inf (PopulationInfo): Founder population
        n_snp (int): No. causal SNPs
        gblup (int): No. neutral SNPs
        n_sample (int): No. individuals

    Retunrs:
        numpy.ndarray: SNP matrix (no. haplotypes x no. SNPs)
    '''
    if gblup is None:
        n_total = n_snp
    else:
        n_total = n_snp + gblup
    rate = (4e-8, 2e-7)  # to simulate almost independent SNPs
    gsa_out = np.empty((0, 2*n_sample), dtype=np.int32)
    gsa_count = 0  # until n_snp
    while True:
        gsa_ts = mp.sim_ancestry(samples=n_sample,
                                 recombination_rate=rate[1],
                                 model=mp.StandardCoalescent(),
                                 sequence_length=200,
                                 population_size=10_000)
        gsa_ts = mp.sim_mutations(gsa_ts,
                                  rate=rate[0],
                                  model=mp.BinaryMutationModel(),
                                  discrete_genome=False,
                                  keep=False)
        for g_ts in gsa_ts.variants():
            tmp_gen = g_ts.genotypes
            # skip monomorphic
            a_fq = np.sum(tmp_gen)
            if a_fq in (0, 2*n_sample):
                continue
            gsa_out = np.append(gsa_out, np.array([tmp_gen]), axis=0)
            gsa_count += 1
            if gsa_count == n_total:
                return gsa_out.T.copy()
# run_msprime


def coalescent_simulation(model, par_inf, n_snp, gblup, n_pop, fst, n_female,
                          n_male):
    '''
    Control the coalescent process

    Args:
        model (str): 'WF' for Wright-Fisher, 'SP' for structed populations
        par_inf (PopulationInfo): Founder population
        n_snp (int): No. causal SNPs
        gblup (int): No. neutral SNPs
        n_pop (int): No. populations
        fst (float): Average Fst value among populations
        n_female (tuple): Nos. females in populations
        n_male (tuple): Nos. males in populations
    '''
    # Wright-Fisher model
    if model == 'WF':
        snp_mat = run_msprime(n_snp, gblup, par_inf.n_f + par_inf.n_m)
    # Structured populations
    elif model == 'SP':
        snp_mat = structured_population(n_snp, gblup, n_pop, fst, n_female,
                                        n_male)
    else:
        sys.exit(f'{model} is under development')
    add_founder_snp(par_inf, snp_mat)
# coalescent_simulation


def main():
    '''
    main
    '''
    print('A module for coalescent simulation')
# main


if __name__ == '__main__':
    main()
