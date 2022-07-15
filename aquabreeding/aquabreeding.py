'''
module for simulating aquabreeding program
version 0.3.0

Copyright (C) 2022 Shohei Takuno

aquabreeding.py allows simulating an aquabreeding program.
First, a founder population with {founder_size} diploid
individuals by coalesent simulation, implemented in msprime.
Assuming all SNPs have some effects on the phenotype,
the progeny population with {progeny_size} individual and
their phenotypes are simulated.
The individuals with large phenotypic values are selected
based on estimated breeding values or phenotypic
values.  The selected ones are used as the founder population
in the next generation.
'''

import sys
import copy
import operator
import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
import msprime


def add_break_point(ab_posi, ab_ls):
    '''
    add recombination break point at the list of a chromatid
    args: recombination break point (int), genome info (list)
    return the index of inserted segment
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
    step 1: get the indices, where crossing-over occurs
    step 2: exchange the segments behind the break point
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


def gametogenesis(pat_chrom, mat_chrom, chromosome_len, cm_mb):
    '''
    a function for simulating gametogenesis
    args: paternal chromosome (list), maternal chromosome (list), chrom len (bp), cM/Mb (float)
    return one of four chromatids (list)
    '''
    # bivalent chromosomes
    pat_1 = [[i_1[0], i_1[1], i_1[2]] for i_1 in pat_chrom]
    pat_2 = [[i_1[0], i_1[1], i_1[2]] for i_1 in pat_chrom]
    mat_1 = [[i_1[0], i_1[1], i_1[2]] for i_1 in mat_chrom]
    mat_2 = [[i_1[0], i_1[1], i_1[2]] for i_1 in mat_chrom]
    # crossing-over
    # expected number of crossing-over events
    # multiplied by 2 to simulate bivalent chromosomes
    ex_n_rec = 2.0*cm_mb*chromosome_len/(10.0**6)/100.0
    if ex_n_rec < 1.0:
        sys.exit('Expected number of crossing-over events is less than 1')
    # '1+' means obligate chiasma
    n_rec = 1 + np.random.poisson(ex_n_rec-1)
    for _ in range(n_rec):
        # recombination break point
        break_point = np.random.randint(1, chromosome_len)+0.5
        # two of which chromatids are recombined
        bivalent_1 = np.random.randint(4)
        if bivalent_1 == 0:
            crossing_over(break_point, pat_1, mat_1)
        elif bivalent_1 == 1:
            crossing_over(break_point, pat_1, mat_2)
        elif bivalent_1 == 2:
            crossing_over(break_point, pat_2, mat_1)
        elif bivalent_1 == 3:
            crossing_over(break_point, pat_2, mat_2)
    bivalent_2 = np.random.randint(4)
    if bivalent_2 == 0:
        return pat_1
    if bivalent_2 == 1:
        return pat_2
    if bivalent_2 == 2:
        return mat_1
    if bivalent_2 == 3:
        return mat_2
    return [-777, -777, -777]
    # gametogenesis


def generate_snp_array(n_sample, n_snp, rate):
    '''
    generate independent SNPs
    return np.array
    rows: haplotype, columns: loci
    '''
    gsa_out = np.empty((0, 2*n_sample), dtype=int)
    gsa_count = 0  # until n_snp
    while True:
        gsa_ts = msprime.sim_ancestry(
            samples=n_sample,
            recombination_rate=rate[1],
            sequence_length=200, population_size=10_000)
        gsa_ts = msprime.sim_mutations(gsa_ts,
            rate=rate[0],
            model=msprime.BinaryMutationModel(),
            discrete_genome=False, keep=False)
        for g_ts in gsa_ts.variants():
            # print (type (g_ts.genotypes))
            tmp_gen = g_ts.genotypes
            # print (tmp_gen)
            gsa_out = np.append(gsa_out, np.array([tmp_gen]), axis=0)
            gsa_count += 1
            if gsa_count == n_snp:
                return gsa_out.T  # np.array
    # generate_snp_array


def two_different_randint(td_max):
    '''
    pickup two different random integer
    (0~td_max-1)
    deprecated later
    '''
    td_i1 = np.random.randint(td_max)
    while True:
        td_i2 = np.random.randint(td_max)
        if td_i1 != td_i2:
            return td_i1, td_i2
    # two_different_randint


class IndividualInfo:
    '''
    class for an individual and its genome
    '''
    def __init__(self, chrom, input_id):
        '''
        storing the information of an individual
        args: tuple(chromosome number, bp, cM/Mb), id
        '''
        self.ind_id = input_id
        self.n_chrom = chrom[0]
        self.l_chrom = chrom[1]
        self.cm_mb = chrom[2]
        # paternal and maternal genomes
        self.chrom_pat = ['' for ii_i in range(self.n_chrom)]
        self.chrom_mat = ['' for ii_i in range(self.n_chrom)]
        for ii_i in range(self.n_chrom):
            # begin segment (bp), end segment (bp), genotype ID in the founder
            self.chrom_pat[ii_i] = [[1, self.l_chrom, 2*self.ind_id]]
            self.chrom_mat[ii_i] = [[1, self.l_chrom, 2*self.ind_id+1]]
        self.father_id = -1
        self.mother_id = -1
        # incorporating sex later
        # 0: male, i: female
        self.sex = -7
        # __init__

    def copy_gametes(self, gamete_1, gamete_2, chrom_id, par_id):
        '''
        exchange chromosome by new gametes
        args: gamete1 (list), gamete2 (list), chromosome ID, tuple(father ID, mother ID)
        '''
        self.chrom_pat[chrom_id] = [[i_1[0], i_1[1], i_1[2]] for i_1 in gamete_1]
        self.chrom_mat[chrom_id] = [[i_1[0], i_1[1], i_1[2]] for i_1 in gamete_2]
        self.father_id = par_id[0]
        self.mother_id = par_id[1]
        # copy_gametes

    def identity_by_descent(self):
        '''
        calculating identity by descent
        '''
        def tmp_ibd(ls_1, ls_2):
            '''
            calculate IBD region in a chromosome
            '''
            # print (ls_1, ls_2)
            pat_i = mat_i = 0
            homo_reg = 0.0
            beg_seg = 1
            end_seg = -777
            while True:
                if ls_1[pat_i][1] < ls_2[mat_i][1]:
                    end_seg = ls_1[pat_i][1]
                    if ls_1[pat_i][2] == ls_2[mat_i][2]:
                        homo_reg += (end_seg - beg_seg + 1)
                    pat_i += 1
                    beg_seg = end_seg + 1
                else:
                    end_seg = ls_2[mat_i][1]
                    if ls_1[pat_i][2] == ls_2[mat_i][2]:
                        homo_reg += (end_seg - beg_seg + 1)
                    mat_i += 1
                    beg_seg = end_seg + 1
                if end_seg == self.l_chrom:
                    break
            if homo_reg > self.l_chrom:
                print(f'===\n{homo_reg}\n{ls_1}\n{ls_2}')
                sys.exit()
            return homo_reg
            # tmp_ibd
        tmp_v = 0.0
        for j in range(self.n_chrom):
            tmp_v += tmp_ibd(self.chrom_pat[j], self.chrom_mat[j])
        tmp_v /= (self.n_chrom * self.l_chrom)
        if tmp_v > 1.0:
            sys.exit('Error in calculating_ibd.  The value is more than 1.')
        return tmp_v
        # ideitiy_by_descent

    def print_individual(self):
        '''
        display-print individual info
        '''
        print(f'Individual ID = {self.ind_id}')
        print(f'  father = {self.father_id}')
        print(f'  mother = {self.mother_id}')
        for pi_j in range(self.n_chrom):
            print(f'    {pi_j}th chromosome')
            print(f'      {self.chrom_pat[pi_j]}')
            print(f'      {self.chrom_mat[pi_j]}')
        # print_individual

    def call_gametogenesis(self, ii_j):
        '''
        call gametogenesis function for {ii_j}th chromosome
        return a recombined chromatid
        '''
        return gametogenesis(self.chrom_pat[ii_j], self.chrom_mat[ii_j], self.l_chrom, self.cm_mb)
        # call_gamatogenesis
    # class IndividualInfo


def progeny_snp_array(snp_array, snp_dict, n_snp, progeny_pop):
    '''
    generating progeny's snp_array
    args: snp_array (founder), snp information, # snp, progeny_pop info
    return np.ndarray
    '''
    def search_genotype(sp_ls, sp_pos):
        '''
        get genotype in a given chromosome position
        '''
        # print (sp_pos)
        # print (sp_ls)
        for sg_ls in sp_ls:
            if sg_ls[0] <= sp_pos <= sg_ls[1]:
                return sg_ls[2]
        return -777
        # search genotype
    n_progeny = len(progeny_pop)
    tmp_snp_array = np.full((2*n_progeny, n_snp), -777777, dtype=int)
    # for each SNP
    for ps_j in range(n_snp):
        ps_chrom = snp_dict[ps_j]['chrom']
        ps_pos = snp_dict[ps_j]['position']
        # for each progeny
        for ps_i in range(n_progeny):
            # get genotype
            ps_pat = search_genotype(progeny_pop[ps_i].chrom_pat[ps_chrom], ps_pos)
            ps_mat = search_genotype(progeny_pop[ps_i].chrom_mat[ps_chrom], ps_pos)
            # print (self.snp_array[ps_pat][ps_j], self.snp_array[ps_pat][ps_j])
            tmp_snp_array[2*ps_i][ps_j] = snp_array[ps_pat][ps_j]
            tmp_snp_array[2*ps_i+1][ps_j] = snp_array[ps_mat][ps_j]
    return tmp_snp_array
    # progeny_snp_array


def convert_snp_to_genotype_array(snp_array):
    '''
    input snp_array, where rows: haplotype, columns: loci
    snp_array is generated by the function generate_snp_array
    SNPs are biallelic
    Genotypes are
    0: 0/0, 1: 0/1, 2: 1/1
    return np.array, where rows: genotypes, columns: loci
    '''
    cs_row, cs_column = np.shape(snp_array)
    cs_out = np.empty((0, cs_column), dtype=int)
    for cs_k in range(0, cs_row, 2):
        cs_tmp = [cs_i+cs_j for cs_i, cs_j in zip(snp_array[cs_k], snp_array[cs_k+1])]
        cs_out = np.append(cs_out, np.array([cs_tmp]), axis=0)
    return cs_out
    # convert_snp_to_genotpe_array


# version 3
def gebv_calculation(vec_y, g_mat, em_threshhold=1e-4):
    '''
    calculating genomic breeding values by solving MME
    Z design matrix in MME is an identity matrix.
    X design matrix conteins only an intercept.
    REML and EM algorithm are applied to estimate variance components
    which follows Saito & Iwaisaki (1995) The EM-REML iteration equations
    for estimation of variance components in a sire and dam model
    Japanses Journal of Biometrics
    args: phenotype vector, genomic/numerator relationship matrix (n x n)
    return beta (p x 1), breeding values (n x 1), additive variance, residual
    speed is very slow due to calculating an inverse matrix and running EM
    algorithm
    '''
    # size of matrix
    # X_rows should be the same as n_sample
    n_sample = len(vec_y)
    var_y = np.var(vec_y, ddof=1)
    vec_y = np.array([vec_y]).T
    des_x = np.array([np.ones(n_sample)]).T
    x_rows, x_columns = np.shape(des_x)
    des_xt = des_x.T
    tmp_var_g = var_y/2.0
    tmp_var_e = var_y/2.0
    tmp_heri = tmp_var_g/(tmp_var_g+tmp_var_e)
    g_inv = np.linalg.inv(g_mat)
    # scale of variance values
    # NOTE: if Z is not an identity matrix, modify here
    var_scale = np.array([[x_rows, 0], [0, x_rows]])
    var_scale = np.linalg.inv(var_scale)
    # var_scale = np.array ([[1.0/X_rows, 0], [0, 1.0/X_rows]])
    tmp_lhs11 = des_xt @ des_x
    tmp_lhs22 = np.identity(x_rows)
    tmp_rhs11 = des_xt @ vec_y
    des_w = np.block([des_x, tmp_lhs22])
    des_wt = des_w.T
    # EM algorithm
    while True:
        # Loop contains two steps
        # 1st: to estimate beta and breeding values given variance components
        # 2nd: to estimate variance components given beta and breeding values
        inv_var_g = 1.0/tmp_var_g
        inv_var_e = 1.0/tmp_var_e
        # 1st step: solve MME
        # left-hand size of MME
        lhs11 = tmp_lhs11 * inv_var_e
        lhs12 = des_xt * inv_var_e
        lhs21 = des_x * inv_var_e
        lhs22 = tmp_lhs22 * inv_var_e + g_inv * inv_var_g
        lhs = np.block([[lhs11, lhs12], [lhs21, lhs22]])
        coef_mat = np.linalg.inv(lhs)
        # right-hand size of MME
        rhs11 = tmp_rhs11 * inv_var_e
        rhs21 = vec_y * inv_var_e
        rhs = np.block([[rhs11], [rhs21]])
        # solve MME
        solution = coef_mat @ rhs
        beta = solution[0: x_columns]
        gebv = solution[x_columns: x_rows+x_columns]
        # 2nd step: variance components estimation
        # following Saito & Iwaisaki (1995)
        error = vec_y - des_x @ beta - gebv
        # C = np.linalg.inv (lhs)
        c_22 = coef_mat[x_columns: x_columns+x_rows, x_columns: x_columns+x_rows]
        tmp_var_e = (error.T @ error + np.trace(des_w @ coef_mat @ des_wt)) * var_scale[1, 1]
        tmp_var_g = (gebv.T @ g_inv @ gebv + np.trace(g_inv @ c_22)) * var_scale[0, 0]
        tmp_var_e = tmp_var_e[0][0]
        tmp_var_g = tmp_var_g[0][0]
        heri = tmp_var_g/(tmp_var_g+tmp_var_e)
        # print (tmp_var_g, tmp_var_e, tmp_var_g/(tmp_var_g+tmp_var_e))
        if abs(tmp_heri - heri) < em_threshhold:
            return beta, gebv, tmp_var_g, tmp_var_e
        tmp_heri = heri
    # gebv_calculation


def numerator_relationshp_matrix(par_dict, progeny_pop, founder_size):
    '''
    calculating numerator relationship matrix
    parental relationship of all founders and current progeny are required
    args: dictionary of parents of founders in all generations
          current progeny populatin info
          size of founder population
    return np.ndarray (n x n)
    '''
    par_key = list(par_dict.keys())
    par_key.sort()
    n_par = len(par_key)
    progeny_size = len(progeny_pop)
    n_mat = np.full((n_par+progeny_size, n_par+progeny_size), -777.0, dtype=float)
    # individuals with unknown parents, i.e., founder
    for i in range(founder_size):
        for j in range(founder_size):
            if i == j:
                n_mat[i][j] = 1.0
            else:
                n_mat[i][j] = 0.0

    # individuals with known parents
    # for j in range (founder_size, n_par):
    for j in range(founder_size, n_mat.shape[1]):
        for i in range(0, j+1):
            pars = pard = -777
            if j < n_par:
                pars = par_dict[j]['pat']
                pard = par_dict[j]['mat']
            else:
                pars = progeny_pop[j-n_par].father_id
                pard = progeny_pop[j-n_par].mother_id
            # print (pars, pard)
            ais = 0
            aid = 0
            if i != j:
                if pars != -1:
                    ais = n_mat[i][pars]
                if pard != -1:
                    aid = n_mat[i][pard]
                n_mat[i][j] = 0.5 * (ais + aid)
                n_mat[j][i] = n_mat[i][j]
            else:
                if pars != -1 and pard != -1:
                    n_mat[i][j] = 1.0 + 0.5 * n_mat[pars][pard]
                else:
                    n_mat[i][j] = 1.0
    return n_mat[n_par:, n_par:]
    # numerator_relationship_matrix


def different_randint_set(r_max, n_set):
    '''
    obtain {n_set} sets of different randint (0 ~ {r_max})
    return list [[s0_1, s0_2], [s1_1m s1_2], [s2_1, s2_2]...]]
    '''
    tmp_rand = list(np.random.randint(r_max, size=(3*n_set, 2)))
    res_rand = [[x[0], x[1]] for x in tmp_rand if x[0] != x[1]]
    len_rand = len(res_rand)
    if len_rand < n_set:
        print(f'{len_rand=} {n_set=}')
        sys.exit('error in different_randint_set')
    return res_rand
    # different_randint_set


class SimAquabreeding:
    '''
    class for aquabreeding program
    '''
    def __init__(self, founder_size, chrom, progeny_size, n_snp,
                       additive_var, residual_var, rate=(2e-8, 4e-7)):
        '''
        generating a founder populations at 0th generation
        storing parental and progeny populations
        args: founder_size: founder population size
              chrom: tuple (chromosome number, bp, cM/Mb)
              progeny_size: population size by parental population mating
              n_snp: the number of snp with additive effects
              rate: tuple (mutation rate, recombination rate)
              additive/residual_var = additive and residual vairance values
              NOTE: additive_var is not acutal additive genetic variance but
                    the variance of effect size of SNPs
        '''
        self.founder_size = founder_size
        self.chromosome = chrom
        # size of population for breeding
        self.progeny_size = progeny_size
        # list for individuals through breeding generations
        self.parental_pop = []
        # list for progeny population
        self.progeny_pop = []
        # snp_info
        self.n_snp = n_snp
        self.rate = rate
        self.var_g = additive_var
        self.var_e = residual_var
        # founder population, id = 0 - self.founder_size-1
        # generate founder
        self.tmp_id = 0
        # dictonary for parents' ID
        self.par_dict = {}
        for _ in range(self.founder_size):
            self.parental_pop.append(IndividualInfo(self.chromosome, self.tmp_id))
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = -1
            self.par_dict[self.tmp_id]['mat'] = -1
            self.tmp_id += 1
        # other variables
        self.snp_array = []
        self.snp_dict = {}
        self.tmp_snp_array = []
        self.effect_size = []
        self.true_breed_value = []
        self.esti_breed_value = []
        self.pheno_value = []
        self.res_beta = []
        self.hat_var_g = -7.0
        self.hat_var_e = -7.0
        # __init__

    def print_population(self, target='founder', n_show=10):
        '''
        display-print a parental/progeny population for debugging
        args: parental_pop ('founder') or progeny_pop ('progeny')
              {n_show} individuals are shown
        '''
        if target == 'founder':
            for pp_i in range(n_show):
                self.parental_pop[pp_i].print_individual()
        elif target == 'progeny':
            for pp_i in range(n_show):
                self.progeny_pop[pp_i].print_individual()
        else:
            sys.exit('target should be founder or progeny')
        # print_population

    def progeny_random_mating(self):
        '''
        generating a progeny population by crossing self.parental_pop
        random mating, ignoreing sex
        progeny's SNP array, true breeding values and  phenotypes
        '''
        # generate breeding population
        self.progeny_pop = []
        randint_set = different_randint_set(self.founder_size, self.progeny_size)
        # generate progenies
        for gp_i in range(self.progeny_size):
            # pickup two individual (father and mother)
            gp_p1 = randint_set[gp_i][0]
            gp_p2 = randint_set[gp_i][1]
            # set progeny's information
            self.progeny_pop.append(IndividualInfo(self.chromosome, -1))
            # meiosis for each chromosome
            for gp_j in range(self.chromosome[0]):
                gamete_1 = self.parental_pop[gp_p1].call_gametogenesis(gp_j)
                gamete_2 = self.parental_pop[gp_p2].call_gametogenesis(gp_j)
                self.progeny_pop[gp_i].copy_gametes(gamete_1, gamete_2, gp_j,
                          (self.parental_pop[gp_p1].ind_id, self.parental_pop[gp_p2].ind_id))
        # generating_progeny_random

    def progeny_partial_diallel_mating(self):
        '''
        generating progeny by partial round-robin mating
        '''
        # messy script with no flexibility... modify this later
        self.progeny_pop = []
        p_count = 0
        parg_1 = int(self.founder_size/2)  # parent / 2
        n_pro_each = int(self.progeny_size / self.founder_size)
        for gp_p1 in range(parg_1):
            for gp_p2 in range(parg_1+gp_p1, parg_1+gp_p1+2):
                if gp_p2 == self.founder_size:
                    gp_p2 = parg_1
                for _ in range(n_pro_each):
                    self.progeny_pop.append(IndividualInfo(self.chromosome, -1))
                    for gp_j in range(self.chromosome[0]):
                        gamete_1 = self.parental_pop[gp_p1].call_gametogenesis(gp_j)
                        gamete_2 = self.parental_pop[gp_p2].call_gametogenesis(gp_j)
                        self.progeny_pop[p_count].copy_gametes(gamete_1, gamete_2, gp_j,
                                                          (self.parental_pop[gp_p1].ind_id,
                                                           self.parental_pop[gp_p2].ind_id))
                    p_count += 1
        if p_count != self.progeny_size:
            print(f'error in progeny_partial_diallel_mating {p_count=} {self.progeny_size=}')
            sys.exit()
        # generating_progeny_partial

    def snp_info(self):
        '''
        generating SNPs in numpy.ndarray for the founder population
        store the all SNP info (allele, effect size, chrom, and position)
        '''
        self.snp_array = generate_snp_array(n_sample=self.founder_size,
                                            n_snp=self.n_snp, rate=self.rate)
        self.effect_size = np.array(np.random.normal(loc=0.0, scale=self.var_g, size=self.n_snp))
        self.snp_dict = {}
        for si_i in range(self.n_snp):
            self.snp_dict[si_i] = {}
            self.snp_dict[si_i]['chrom'] = np.random.randint(self.chromosome[0])
            self.snp_dict[si_i]['position'] = np.random.randint(self.chromosome[1]) + 1
        # snp_info

    def breeding_value_nrm(self, g_blup=True):
        '''
        calculating breeding values with numerator relationship matrix
        together wih generatin snp array, genotype array, phenotypes
        if g_blup is false, numerator relationship matrix, breeding values,
        variance components are not calculated
        '''
        self.tmp_snp_array = progeny_snp_array(self.snp_array, self.snp_dict, self.n_snp,
                                               self.progeny_pop)
        genotype_array = convert_snp_to_genotype_array(self.tmp_snp_array)
        # calculating true breeding value and phenotypic value
        self.true_breed_value = np.empty(self.progeny_size)
        self.pheno_value = np.empty(self.progeny_size)
        rand_norm = np.random.normal(loc=0.0, scale=self.var_e, size=self.progeny_size)
        for bvg_i in range(self.progeny_size):
            self.true_breed_value[bvg_i] = np.dot(genotype_array[bvg_i], self.effect_size.T)
            self.pheno_value[bvg_i] = self.true_breed_value[bvg_i] + rand_norm[bvg_i]
        # numerator relationsip matrix and genomic breeding value
        if g_blup:
            n_mat = numerator_relationshp_matrix(self.par_dict, self.progeny_pop, self.founder_size)
            self.res_beta, self.esti_breed_value, self.hat_var_g, self.hat_var_e = gebv_calculation(
                                                  self.pheno_value, n_mat, em_threshhold=1e-4)
        # breeding_value_nrm

    def select_by_breeding_value(self):
        '''
        generate next parental population based on estimated breeding values
        '''
        id_bv = [[i, self.esti_breed_value[i]] for i in range(self.progeny_size)]
        id_bv.sort(reverse=True, key=operator.itemgetter(1))
        for i in range(self.founder_size):
            self.parental_pop[i] = copy.deepcopy(self.progeny_pop[id_bv[i][0]])
            self.parental_pop[i].ind_id = self.tmp_id
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = self.parental_pop[i].father_id
            self.par_dict[self.tmp_id]['mat'] = self.parental_pop[i].mother_id
            self.tmp_id += 1
        # select_by_breeding_value

    def select_by_pedigree(self):
        '''
        generate next parental population based on pedigree
        This function is used together with progeny_partial_diallel_mating
        '''
        id_bv = []
        for i in range(self.progeny_size):
            pat_i = self.progeny_pop[i].father_id
            mat_i = self.progeny_pop[i].mother_id
            if pat_i < mat_i:
                id_bv.append([i, pat_i, mat_i, self.esti_breed_value[i]])
            else:
                id_bv.append([i, mat_i, pat_i, self.esti_breed_value[i]])
        id_bv.sort(reverse=True, key=operator.itemgetter(3))
        # search for the individual with the largest breeding value
        # within subpopulation (parents are the same)
        # ugly... modify this later
        for i in range(self.progeny_size-1):
            if id_bv[i][0] == -7:
                continue
            for j in range(i+1, self.progeny_size):
                if id_bv[j][0] == -7:
                    continue
                if id_bv[i][1] == id_bv[j][1] and id_bv[i][2] == id_bv[j][2]:
                    id_bv[j][0] = -7
        # print (id_bv)
        p_count = 0
        for id_list in id_bv:
            if id_list[0] == -7:
                continue
            self.parental_pop[p_count] = copy.deepcopy(self.progeny_pop[id_list[0]])
            self.parental_pop[p_count].ind_id = self.tmp_id
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = self.parental_pop[p_count].father_id
            self.par_dict[self.tmp_id]['mat'] = self.parental_pop[p_count].mother_id
            self.tmp_id += 1
            p_count += 1
        if p_count != self.founder_size:
            print(f'the number of new parental pop is wrong {p_count}')
            sys.exit()
        # select_by_pedigree

    def select_at_random(self):
        '''
        no selection
        no need to calculate breeding values
        '''
        # no selection
        rand_rand = np.random.rand(self.progeny_size)
        id_bv = [[i, rand_rand[i]] for i in range(self.progeny_size)]
        id_bv.sort(reverse=True, key=operator.itemgetter(1))
        for i in range(self.founder_size):
            self.parental_pop[i] = copy.deepcopy(self.progeny_pop[id_bv[i][0]])
            self.parental_pop[i].ind_id = self.tmp_id
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = self.parental_pop[i].father_id
            self.par_dict[self.tmp_id]['mat'] = self.parental_pop[i].mother_id
            self.tmp_id += 1
        # select_at_random

    def select_by_phenotype(self):
        '''
        select based only on phenotypes
        no need to calculate breeding values
        '''
        # selection
        id_bv = [[i, self.pheno_value[i]] for i in range(self.progeny_size)]
        id_bv.sort(reverse=True, key=operator.itemgetter(1))
        for i in range(self.founder_size):
            self.parental_pop[i] = copy.deepcopy(self.progeny_pop[id_bv[i][0]])
            self.parental_pop[i].ind_id = self.tmp_id
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = self.parental_pop[i].father_id
            self.par_dict[self.tmp_id]['mat'] = self.parental_pop[i].mother_id
            self.tmp_id += 1
        # select_by_phenotype

    def calculate_ibd(self):
        '''
        calculating the probability of IBD
        '''
        ibd_res = []
        for i in range(self.progeny_size):
            tmp_v2 = self.progeny_pop[i].identity_by_descent()
            ibd_res.append(tmp_v2)
        return ibd_res
    # class SimAquabreeding


def main():
    '''
    main
    '''
    print('module for simulating aquabreeding progrem')


if __name__ == '__main__':
    main()
















