'''
module for simulating aquabreeding program
version 0.5.0

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

from aquabreeding import popgen as pg
from aquabreeding import gametogenesis as gg
from aquabreeding import gblup as gb


def different_randint_set(r_max, n_set):
    '''
    obtain {n_set} sets of different randint (0 ~ {r_max})
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


class ChromInfo:
    '''
    class for a chromosome
    args: chrom: tuple(chrom num, chrom len, cM/Mb)
    input_id: ID of parental individuals
    '''
    def __init__(self, chrom, input_id):
        self.chrom = chrom
        self.c_id = input_id
        # begin segment (bp), end segment (bp), genotype ID in the founder
        self.ch_pat = np.array([[1, self.chrom[1], 2*self.c_id]])
        self.ch_mat = np.array([[1, self.chrom[1], 2*self.c_id+1]])
        # __init__

    def call_gametogenesis2(self, ex_n_rec):
        '''
        call gametogenesis in gametogenesis.py
        '''
        return gg.gameto_genesis(self.ch_pat,
                                 self.ch_mat,
                                 ex_n_rec,
                                 self.chrom[1])
# class ChromInfo


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
        self.chrom = chrom
        # paternal and maternal genomes
        self.ch_ls = []
        for _ in range(self.chrom[0]):
            self.ch_ls.append(ChromInfo(self.chrom, self.ind_id))
        # ID of parents
        self.pat_id = -1
        self.mat_id = -1
        # incorporating sex later
        # 0: male, i: female
        self.sex = -7
        # __init__

    def print_individual(self):
        '''
        display-print individual info
        '''
        print('======')
        print(f'Individual ID = {self.ind_id}')
        print(f'  father = {self.pat_id}')
        print(f'  mother = {self.mat_id}')
        for pi_j in range(self.chrom[0]):
            print(f'{pi_j}th chromosome')
            print('paternal')
            print(f'{self.ch_ls[pi_j].ch_pat}')
            print('maternal')
            print(f'{self.ch_ls[pi_j].ch_mat}')
        # print_individual

    def call_gametogenesis(self, c_j, ex_n_rec):
        '''
        call gametogenesis function for {ii_j}th chromosome
        return a recombined chromatid
        '''
        return self.ch_ls[c_j].call_gametogenesis2(ex_n_rec)
        # call_gamatogenesis

    def copy_gametes(self, gamete_1, gamete_2, ch_id):
        '''
        exchange chromosome by new gametes
        args: gamete1 (list), gamete2 (list), chromosome ID,
        '''
        self.ch_ls[ch_id].ch_pat = gamete_1.copy()
        self.ch_ls[ch_id].ch_mat = gamete_2.copy()
        # copy_gametes
# class IndividualInfo


class SimAquabreeding:
    '''
    class for aquabreeding program
    '''
    def __init__(self, founder_size, chrom, progeny_size, n_snp,
                 additive_var, residual_var):
        '''
        generating a founder populations
        storing parental and progeny populations
        args: founder_size: founder population size
              chrom: tuple(chromosome number, bp, cM/Mb)
              progeny_size: population size by parental population mating
              n_snp: the number of snp with additive effects
              additive/residual_var = additive and residual vairance values
              NOTE: additive_var is not acutal additive genetic variance but
                    the variance of effect size of SNPs
        '''
        # no. parents
        self.n_founder = founder_size
        self.chrom = chrom
        # no. progeny
        self.n_progeny = progeny_size
        # Population used as parents
        self.par_pop = []
        # Progeny population
        self.pro_pop = []
        # snp_info
        self.n_snp = n_snp
        self.v_g = additive_var
        self.v_e = residual_var
        # founder population, id = 0 - self.founder_size-1
        # generate founder
        self.tmp_id = 0  # ID of founders through generations
        # dictonary for parents' ID
        self.par_dict = {}
        for _ in range(self.n_founder):
            self.par_pop.append(IndividualInfo(self.chrom,
                                               self.tmp_id))
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = -1
            self.par_dict[self.tmp_id]['mat'] = -1
            self.tmp_id += 1
        # other variables
        self.par_snp = []  # SNP data of the founder pop
        self.snp_dict = {}  # SNP info (chrom and position)
        self.pro_snp = []  # SNP data of the progeny pop
        self.effect_size = []  # effect size of SNPs
        self.true_bv = []  # true breeding value
        self.hat_bv = []   # predicted breeding value
        self.pheno_v = []  # phenotypic value
        self.hat_beta = []  # estimated beta in ME
        self.hat_vg = -7.0  # esimated additive variance
        self.hat_ve = -7.0  # estimated residual variance
        # the expected number of crossing over evtens
        # this value must be equall to or bigger than 1
        self.ex_n_rec = 2.0*self.chrom[2]*self.chrom[1]/(10.0**6)/100.0
        if self.ex_n_rec < 1.0:
            sys.exit('Expected number of crossing-over events is less than 1')
        # __init__

    def print_population(self, target='founder', n_show=10):
        '''
        display-print a parental/progeny population for debugging
        args: target='founder' for parental population
                      or
                     'progeny' for progeny population
              n_show no. of displayed individuals
        '''
        if target == 'founder':
            if n_show > self.n_founder:
                print('Error in print_population')
                sys.exit('n_show is bigger than the founder size')
            for pp_i in range(n_show):
                self.par_pop[pp_i].print_individual()
        elif target == 'progeny':
            if n_show > self.n_progeny:
                print('Error in print_population')
                sys.exit('n_show is bigger than the progeny size')
            for pp_i in range(n_show):
                self.pro_pop[pp_i].print_individual()
        else:
            sys.exit('target should be founder or progeny')
        # print_population

    def snp_info(self, rate=(2e-8, 4e-7)):
        '''
        generating SNPs in ny.ndarray for the founder population
        store the all SNP info (allele, effect size, chrom, and position)
        assuming panmictic population in this version
        ust scripts in popgen.py
        '''
        self.par_snp = pg.generate_snp(n_sample=self.n_founder,
                                       n_snp=self.n_snp, rate=rate)
        self.effect_size = np.array(np.random.normal(loc=0.0,
                                                     scale=self.v_g,
                                                     size=self.n_snp))
        self.snp_dict = {}
        rand_c = np.random.randint(self.chrom[0], size=self.n_snp)
        rand_p = np.random.randint(self.chrom[1], size=self.n_snp) + 1
        for si_i in range(self.n_snp):
            self.snp_dict[si_i] = {}
            self.snp_dict[si_i]['chrom'] = rand_c[si_i]
            self.snp_dict[si_i]['pos'] = rand_p[si_i]
        # snp_info

    def random_mating(self):
        '''
        generating a progeny population by crossing self.par_pop
        random mating, ignoring sex
        progeny's SNP array, true breeding values and phenotypes
        using gametogenesis.py
        '''
        # generate breeding population
        self.pro_pop = []
        randint_set = different_randint_set(self.n_founder,
                                            self.n_progeny)
        # generate progenies
        for g_i in range(self.n_progeny):
            # pickup two individual (father and mother)
            g_p1 = randint_set[g_i][0]
            g_p2 = randint_set[g_i][1]
            # set progeny's information
            self.pro_pop.append(IndividualInfo(self.chrom, -1))
            # gametogenesis for each chromosome
            for g_j in range(self.chrom[0]):
                gamete_1 = self.par_pop[g_p1].call_gametogenesis(g_j,
                                                                 self.ex_n_rec)
                gamete_2 = self.par_pop[g_p2].call_gametogenesis(g_j,
                                                                 self.ex_n_rec)
                self.pro_pop[g_i].copy_gametes(gamete_1,
                                               gamete_2,
                                               g_j)
                self.pro_pop[g_i].pat_id = self.par_pop[g_p1].ind_id
                self.pro_pop[g_i].mat_id = self.par_pop[g_p2].ind_id
        # random_mating

    def breeding_value_nrm(self, g_blup=True, use_jit=True):
        '''
        calculating breeding values with numerator relationship matrix
        together wih generating snp array, genotype array, phenotypes
        if g_blup is false, numerator relationship matrix, breeding values,
        variance components are not calculated
        if use_jit is true, a function with numba jit used to calculate
                            numerator relationship matrix.
                      false, a function written in C++ is used.
        using gblup.py
        '''
        # genotypes of the progeny population
        self.pro_snp = pg.progeny_snp(self.par_snp, self.snp_dict, self.n_snp,
                                      self.pro_pop, self.n_progeny)
        genotype_array = pg.genotype_array(self.pro_snp)
        # calculating true breeding value and phenotypic value
        rand_norm = np.random.normal(loc=0.0,
                                     scale=self.v_e,
                                     size=self.n_progeny)
        self.true_bv = np.dot(genotype_array, self.effect_size.T)
        self.pheno_v = self.true_bv + rand_norm
        # numerator relationsip matrix and genomic breeding value
        if g_blup:
            if use_jit:
                n_mat = gb.nrm_jit(self.par_dict, self.pro_pop, self.n_founder)
            else:
                n_mat = gb.nrm_cpp(self.par_dict, self.pro_pop, self.n_founder)
            (self.hat_beta,
             self.hat_bv,
             self.hat_vg,
             self.hat_ve) = gb.gebv_calculation(self.pheno_v, n_mat)
            # print (self.hat_beta, np.mean (self.pheno_v),
            #        np.mean(self.true_bv))
        else:
            pass
        # breeding_value_nrm

    def select_by_breeding_value(self):
        '''
        generate next parental population based on predicted breeding values
        '''
        id_bv = [[i, self.hat_bv[i]] for i in range(self.n_progeny)]
        id_bv.sort(reverse=True, key=operator.itemgetter(1))
        for i in range(self.n_founder):
            self.par_pop[i] = copy.deepcopy(self.pro_pop[id_bv[i][0]])
            self.par_pop[i].ind_id = self.tmp_id
            self.par_dict[self.tmp_id] = {}
            self.par_dict[self.tmp_id]['pat'] = self.par_pop[i].pat_id
            self.par_dict[self.tmp_id]['mat'] = self.par_pop[i].mat_id
            self.tmp_id += 1
        # select_by_breeding_value

# class SimAquabreeding


def main():
    '''
    main
    '''
    print('The module for simulating aquatic breeding')
    # main


if __name__ == '__main__':
    main()
