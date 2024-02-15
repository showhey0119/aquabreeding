'''
Module for simulating aquaculture breeding

Aquabreeding allows simulating aquaculture breeding.

A founder population is generated together with coalesent simulation,
implemented in msprime.

Progenies are produced following some implemented mating schemes or
user's own setting.

The individuals with large phenotypic/breeding values are selected.
Mass selection, within-family selection, and family selection are implemented.
The selected ones are used as parentsin the next generation.

Phenotype, true/estimated breeding value, inbreeding coefficient,
and variance components are output.

Note:
    This is a beta version.

Todo:
    * Dominance and epistasis
    * Multiple phenotypes
    * Founders from a population with exponential growth or bottleneck
'''

import sys
import numpy as np
from aquabreeding import coalescent as co
from aquabreeding import mating as mt
from aquabreeding import phenotype as ph
from aquabreeding import popinfo as po
from aquabreeding import selection as se


def check_tuple(obj, tag, l_tuple):
    '''
    Check the arguments of AquaBreeding class

    Check if an argument is tuple and
    if the length of the tuple is correct

    Args:
        obj (unknown): A focal variable
        tag (str): Name of the variable
        l_tuple (int): Correct tuple length
    '''
    if not isinstance(obj, tuple):
        sys.exit(f'{tag} should be tuple')
    if len(obj) != l_tuple:
        if tag in ('founder_size',  'progeny_size'):
            e_mess = '(No. males, No. females)'
        if tag == 'chrom':
            e_mess = '(Chrom no., chrom len (bp), male cM/Mb, female cM/Mb)'
        if tag in ('n_female', 'n_male'):
            e_mess = 'Length should be equal to be n_pop'
        sys.exit(f'Length of {tag} is incorrect\n{e_mess}')
# check_tuple


class AquaBreeding:
    '''
    Class for aquaculture breeding

    Args:
        founder_size (tuple): Nos. females and  males in the founder population
        chrom (tuple): Chrom num, chrom len (bp), female cM/Mb, male cM/Mb
        mean_phenotype (float): Mean phenotype
        var_phenotype (float): Variance of phenotype
        heritability (float): Heeritability
        n_wild (int): No. individuals in a extra natural population,
                     default None

    Attributes:
        chrom (tuple): Chrom num, chrom len (bp), female cM/Mb, male cM/Mb
        par_inf (PopulationInfo): Founder/parental population
        pro_inf (PopulationInfo): Progeny population
        n_snp (int): No. causal SNPs
        phe_inf (PhenotypeInfo): Phenotype information
        cross_inf (numpy.ndarray): Index pairs of female and male parents
        gblup (int): No. neutral SNPs
    '''
    def __init__(self, founder_size, chrom, mean_phenotype, var_phenotype,
                 heritability, n_wild=None):
        '''
        constructor
        '''
        # check argument
        check_tuple(founder_size, 'founder_size', 2)
        check_tuple(chrom, 'chrom', 4)
        # chromosome info
        self.chrom = chrom
        # parental population
        self.par_inf = po.PopulationInfo(founder_size, self.chrom, n_wild)
        self.par_inf.init_founder()
        # progeny population
        self.pro_inf = None
        # Phenotype info
        self.phe_inf = ph.PhenotypeInfo(mean_phenotype, var_phenotype,
                                        heritability)
        # No. causal SNPs
        self.n_snp = None
        # No. neutral SNPs
        self.gblup = None
        # Mating design
        self.cross_inf = None
    # __init__

    def snp(self, model, n_snp, gblup=None, n_pop=None, fst=None,
            n_female=None, n_male=None):
        '''
        Generate SNPs

        Args:
            model (str): 'WF' or 'SP'
            n_snp (int): No. causal SNPs
            gblup (int): No. neutral SNPs
            n_pop (int): No. populations
            fst (float): Average Fst value among populations
            n_female (tuple): Nos. females in populations
            n_male (tuple): Nos. males in populations
        '''
        # check args
        if n_pop is not None:
            check_tuple(n_female, 'n_female', n_pop)
            check_tuple(n_male, 'n_male', n_pop)
            if sum(n_female) != self.par_inf.n_f:
                sys.exit(f'Sum of n_female is not {self.par_inf.n_f}')
            if sum(n_male) != self.par_inf.n_m:
                sys.exit(f'Sum of n_male is not {self.par_inf.n_m}')
        self.n_snp = n_snp
        self.gblup = gblup
        # coalescent simulation
        co.coalescent_simulation(model, self.par_inf, self.n_snp, self.gblup,
                                 n_pop, fst, n_female, n_male)
    # snp

    def mating_design(self, design, select_size=None):
        '''
        Set mating design

        Args:
            design (unknown): If str, design should be '1x(int)' for partial
                              factorial cross or 'full' factorial mating.
                              If numpy.ndarray, design contains two columns:
                              index of female parents and index of male
                              parents.
            select_size (tuple): Set nos. selected females and males
        '''
        if select_size is None:
            select_size = (self.par_inf.n_f, self.par_inf.n_m)
        else:
            check_tuple(select_size, 'select_size', 2)
        self.cross_inf = mt.set_mating_design(design, select_size)
    # mating_design

    def mating(self, progeny_size):
        '''
        Mate founder/parental individuals to produce progenies

        Args:
            progeny_size (tuple): Nos. female and male progenies
        '''
        check_tuple(progeny_size, 'progeny_size', 2)
        if self.pro_inf is None:
            self.pro_inf = po.PopulationInfo(progeny_size, self.chrom, None)
            self.pro_inf.init_progeny()
        else:
            self.pro_inf.change_size(progeny_size)
        mt.start_mating(self.cross_inf, self.par_inf, self.pro_inf)
    # mating

    def breeding_value(self, target='BLUP'):
        '''
        Calculate phenotype and breeding value

        Args:
            target (str): If 'BLUP', numerator relationship matrix is used
                          to estimate breeding values.  If 'GBLUP', genomic
                          relationship matrix is used.  If 'no', breeding
                          values are not estimated. Default 'BLUP'
        '''
        # genotype matrix
        self.pro_inf.genotype_matrix(self.n_snp, self.gblup)
        # Calculate phenotype and breeding value
        self.phe_inf.calculate_bv(target, self.par_inf, self.pro_inf,
                                  self.n_snp, self.gblup)
    # breeding_value

    def selection(self, target='bv', method='mass', top_prop=1.0,
                  n_family=-1, select_size=None, max_f=None):
        '''
        Select parents of next generation

        Args:
            target (str): Selection based on breeding value ('bv'),
                          phenotypic value ('phenotype'), or random
                          ('random')
            method (str): How to select from progenies such as mass
                          selection ('mass'), within-family selection
                          ('within-family'),  family selection ('family'),
                          or selection based on A matrix ('FvalueA') or
                          G matrix ('FvalueG')
            top_prop (float): Select progenies with top X% breeding values
                              in within-family selection. Set 0.0 < top_prop
                              <= 1.0.
            n_family (int): No. families to be selected
            select_size (tulple): Number of selected founders, default: None
            max_f (float): F values among selected individuals are less than
                           max_f
        '''
        if select_size is not None:
            check_tuple(select_size, 'select_size', 2)
            self.par_inf.change_size(select_size)
        else:
            select_size = (self.par_inf.n_f, self.par_inf.n_m)
        se.start_selection(self.par_inf, self.pro_inf, self.phe_inf,
                           target, method, self.cross_inf, top_prop,
                           n_family, select_size, max_f)
    # selection

    def get_ibd(self):
        '''
        Output inbreeding coefficient

        Returns:
            numpy.ndarray: Inbreeding coefficient
        '''
        return np.diag(self.pro_inf.a_mat) - 1.0
    # get_ibd

    def get_phenotype(self):
        '''
        Output phenotype

        Returns:
            numpy.ndarray: Phenotype
        '''
        return self.phe_inf.pheno_v
    # get_phenotype

    def get_true_bv(self):
        '''
        Output true breeding value

        Returns:
            numpy.ndarray: True breeding value
        '''
        return self.phe_inf.true_bv
    # get_true_bv

    def get_ebv(self):
        '''
        Output estimated breeding value

        Returns:
            numpy.ndarray: Estimated breeding value
        '''
        return self.phe_inf.hat_bv
    # get_ebv

    def variance_component(self):
        '''
        Outut estimated additive and residual variance

        Returns:
            tuple: Additive and residual variance
        '''
        return self.phe_inf.hat_vg, self.phe_inf.hat_ve
    # variance_component

    def wild_parent(self, num):
        '''
        Replance parents with wild individuals before mating

        Args:
            num (tuple): No. famale/male parents relpaced by wild individuals
        '''
        self.par_inf.replace_wild(num)
    # use_natural
# AquaBreeding


def main():
    '''
    main
    '''
    print('A module for simulating aquacuture breeding')
# main


if __name__ == '__main__':
    main()
