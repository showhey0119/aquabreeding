'''
module for simulating aquaculture breeding

version 0.7.0

aquabreeding allows simulating an aquaculture breeding program.

A founder population is generated together with coalesent simulation,
implemented in msprime.

Progenies are generated following some implemented mating scheme or
user's own setting.

The individuals with large phenotypic/breeding values are selected.
Mass selection or within-family selection are implemented.
The selected ones are used as the parental population in the next
generation.

Phenotype, true/estimated breeding value, inbreeding coefficient,
and variance components are output.

Todo:
* Implement GBLUP
* Implement multiple phenotypes
* Implement founders from structured populations
'''

import sys
import bisect
import numpy as np
from aquabreeding import popgen as pg
from aquabreeding import mating as mt
from aquabreeding import gametogenesis as gg
from aquabreeding import blup as bl
from aquabreeding import selection as se


class ChromInfo:
    '''
    Class for one of the chromosomes of an individual

    Args:
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
        input_id (int): ID (only used for a founder)

    Attributes:
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
        chrom_mat (list): End points of rec segments in maternal chromosome
        chrom_pat (list): End points of rec segments in paternal chromosome
        genotype_mat (list): genotype of maternal chromosome in a founder
        genotype_pat (list): genotype of paternal chromosome in a founder
    '''
    def __init__(self, chrom, input_id):
        '''
        constractor
        '''
        self.chrom = chrom
        # End point of recombination segments
        self.chrom_mat = [self.chrom[1]]
        self.chrom_pat = [self.chrom[1]]
        # genotype of recombination segments
        self.genotype_mat = [2*input_id]
        self.genotype_pat = [2*input_id+1]
    # __init__

    def get_genotype_id(self, snp_pos):
        '''
        Get genotypes in rec segment with SNP

        Args:
            snp_pos (int): SNP position

        Returns:
            genotypes

            - int: genotype in maternal chromosome
            - int: genotype in paternal chromosome
        '''
        index_m = bisect.bisect_left(self.chrom_mat, snp_pos)
        gen_mat = self.genotype_mat[index_m]
        index_p = bisect.bisect_left(self.chrom_pat, snp_pos)
        gen_pat = self.genotype_pat[index_p]
        return gen_mat, gen_pat
    # get_genotype_id
# class ChromInfo


class IndividualInfo:
    '''
    class for an individual and its genome

    Args:
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
        input_id (int): ID of an individual
        pat_id (int): ID of father, default = -1
        mat_id (int): ID of mother, defautl = -1

    Attributes:
        ind_id (int): ID of the individual (only used for parental population)
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
        chrom_ls (list): list of ChromInfo class
        mat_id (int): ID of mother
        pat_id (int): ID of father
        n_recf (float): expected num crossing-over event in female
        n_recm (float): expected num crossing-over event in male
    '''
    def __init__(self, chrom, input_id, pat_id=-1, mat_id=-1):
        '''
        constructor
        '''
        self.ind_id = input_id
        self.chrom = chrom
        # paternal and maternal genomes
        self.chrom_ls = []
        for _ in range(self.chrom[0]):
            self.chrom_ls.append(ChromInfo(self.chrom, self.ind_id))
        # ID of parents (-1 means unknown)
        self.pat_id = pat_id
        self.mat_id = mat_id
        # expected num crossing-over event
        self.n_recf = 2.0*self.chrom[2]*self.chrom[1]/(10.0**6)/100.0
        self.n_recm = 2.0*self.chrom[3]*self.chrom[1]/(10.0**6)/100.0
    # __init__

    def print_individual(self):
        '''
        display-print individual info
        '''
        print('======')
        print(f'Individual ID = {self.ind_id}')
        print(f'  mother = {self.mat_id}')
        print(f'  father = {self.pat_id}')
        for i in range(self.chrom[0]):
            print(f'  {i}th chromosome')
            print('    maternal')
            for j, e_p in enumerate(self.chrom_ls[i].chrom_mat):
                g_p = self.chrom_ls[i].genotype_mat[j]
                if j == 0:
                    print(f'      [1 {e_p:,} {g_p}]')
                else:
                    e_p2 = self.chrom_ls[i].chrom_mat[j-1] + 1
                    print(f'      [{e_p2:,} {e_p:,} {g_p}]')
            print('    paternal')
            for j, e_p in enumerate(self.chrom_ls[i].chrom_pat):
                g_p = self.chrom_ls[i].genotype_pat[j]
                if j == 0:
                    print(f'      [1 {e_p:,} {g_p}]')
                else:
                    e_p2 = self.chrom_ls[i].chrom_pat[j-1] + 1
                    print(f'      [{e_p2:,} {e_p:,} {g_p}]')
    # print_individual

    def get_ibd2(self):
        '''
        Calulate inbreeding coefficient in an individual
        '''
        tmp_v = 0.0
        for j in range(self.chrom[0]):
            tmp_v += pg.identity_by_descent(self.chrom_ls[j])
        tmp_v /= (self.chrom[0] * self.chrom[1])
        if tmp_v > 1.0:
            sys.exit('The value is more than 1 in get_ibd2')
        return tmp_v
    # get_ibd2

    def gamete(self, chrom_id, tag):
        '''
        Produce gamete

        call gametogenesis.py

        Args:
            chrom_id (int): chromosome ID
            tag (int): 0 for female, 1 for male

        Returns:
            new gamete

            - list: recombined chromatid,
            - list: genotypes
        '''
        if tag == 0:  # female
            n_rec = self.n_recf
        elif tag == 1:  # male
            n_rec = self.n_recm
        else:
            sys.exit('tag in gamete should be 0 or 1')
        new_gamete, new_geno = gg.gameto_genesis(self.chrom_ls[chrom_id],
                                                 n_rec)
        return new_gamete, new_geno
    # gamate

    def copy_gametes(self, gamete_1, gamete_2, geno_1, geno_2,
                     chrom_id):
        '''
        exchange chromosome by new gametes

        Args:
            gamete_1 (list): new chromatid from female
            gamete_2 (list): new chromatid from male
            chrom_id (int): chromosome ID

        Note:
            Check whether call-by-reference do something wrong
        '''
        self.chrom_ls[chrom_id].chrom_mat = gamete_1
        self.chrom_ls[chrom_id].chrom_pat = gamete_2
        self.chrom_ls[chrom_id].genotype_mat = geno_1
        self.chrom_ls[chrom_id].genotype_pat = geno_2
    # copy_gametes

    def get_genotype(self, chrom_id, snp_pos):
        '''
        Get genotpye of rec segment that has the SNP

        Args:
            chrom_id (int): chromosome ID
            snp_pos: (int): Position of SNP

        Returns:
            genotypes

            - int: Genotype in maternal chrom
            - int: Genotype in paternal chrom
        '''
        # maternal chrom
        gen_mat, gen_pat = self.chrom_ls[chrom_id].get_genotype_id(snp_pos)
        return gen_mat, gen_pat
    # get_genotype
# class IndividualInfo


class PopulationInfo:
    '''
    class for breeding population information

    Args:
        pop_size (tuple): no. females, no. males (int, int)
                          in a population
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
                       (int, int, float, float)

    Attributes:
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
        n_popf (int): No. of females in a population
        n_popm (int): No. of males in a population
        pop_f (list): list of female IndividualInfo class
        pop_m (list): list of male IndividualInfo class
        tmp_id (int): ID of all parents (only used for a founder)
        d_par (dict): Dict of parents of founders (only used for a founder)
    '''
    def __init__(self, pop_size, chrom):
        '''
        constructur
        '''
        self.chrom = chrom
        self.n_popf = pop_size[0]
        self.n_popm = pop_size[1]
        self.pop_f = []
        self.pop_m = []
        self.tmp_id = 0
        self.d_par = {}
    # __init__

    def init_founder(self):
        '''
        Initiate founder/parental population
        '''
        # female
        for _ in range(self.n_popf):
            self.pop_f.append(IndividualInfo(self.chrom, self.tmp_id))
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = -1
            self.d_par[self.tmp_id]['pat'] = -1
            self.tmp_id += 1
        # male
        for _ in range(self.n_popm):
            self.pop_m.append(IndividualInfo(self.chrom, self.tmp_id))
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = -1
            self.d_par[self.tmp_id]['pat'] = -1
            self.tmp_id += 1
    # init_founder

    def new_founder_id(self):
        '''
        Set IDs for new founder
        '''
        # female
        for i in range(self.n_popf):
            self.pop_f[i].ind_id = self.tmp_id
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = self.pop_f[i].mat_id
            self.d_par[self.tmp_id]['pat'] = self.pop_f[i].pat_id
            self.tmp_id += 1
        # male
        for i in range(self.n_popm):
            self.pop_m[i].ind_id = self.tmp_id
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = self.pop_m[i].mat_id
            self.d_par[self.tmp_id]['pat'] = self.pop_m[i].pat_id
            self.tmp_id += 1
    # new_founder_id

    def init_progeny(self):
        '''
        Initiate progeny population
        '''
        # female
        for _ in range(self.n_popf):
            self.pop_f.append(IndividualInfo(self.chrom, -1))
        # male
        for _ in range(self.n_popm):
            self.pop_m.append(IndividualInfo(self.chrom, -1))
    # init_progeny

    def print_popinfo(self, n_show):
        '''
        Display-print population info

        Args:
            n_show (int): no. displayed individuals (female/male each)
        '''
        if n_show > self.n_popm or n_show > self.n_popf:
            sys.exit('n_show is larger than the no. of female/male')
        print('********** female **********')
        for i in range(n_show):
            self.pop_f[i].print_individual()
        print('\n********** male **********')
        for i in range(n_show):
            self.pop_m[i].print_individual()
    # print_self
# class PopulationInfo


class SNPInfo:
    '''
    class for SNP information

    Args:
        n_snp (int): No. SNPs with effect size
        effect_var (float): Variance of effect size
        founder_size (tuple): Nos. females/males in a founder (int)
        progeny_size (tuple): Nos. females/males in a progeny (int)
        chrom (tuple): chrom num, chrom len, male cM/Mb, female cM/Mb

    Attributes:
        n_snp (int): No. SNPs with effect size
        n_founder (int): No. founders
        par_snp (ndarray): SNP matrix in a founder population
                           (rows: haplotype, column: loci)
        pro_snp (ndarray): SNP matrix in a progeny population
                           (rows: haplotype, column: loci)
        gen_mat (ndarray): genotype matrix
                           (rows: genotype, column: loci)
        effect_size (ndarray): effect size of SNPs
        snp_dict (dict): chrom and position of SNPs
    '''
    def __init__(self, n_snp, effect_var, founder_size, progeny_size,
                 chrom):
        '''
        constructor
        '''
        self.n_snp = n_snp
        self.n_founder = founder_size[0] + founder_size[1]
        self.par_snp = []
        pro_row = progeny_size[0] + progeny_size[1]
        self.gen_mat = np.full((pro_row, n_snp), -7, dtype=np.float64)
        pro_row = 2 * (progeny_size[0] + progeny_size[1])
        self.pro_snp = np.full((pro_row, n_snp), -7, dtype=np.int64)
        if effect_var != -1:
            self.effect_size = np.random.normal(loc=0.0,
                                                scale=np.sqrt(effect_var),
                                                size=self.n_snp)
        self.snp_dict = {}
        # which chromosome
        rand_c = np.random.randint(chrom[0], size=self.n_snp)
        # where within a chromosome
        rand_p = np.random.randint(chrom[1], size=self.n_snp) + 1
        for i in range(self.n_snp):
            self.snp_dict[i] = {}
            self.snp_dict[i]['chrom'] = rand_c[i]
            self.snp_dict[i]['pos'] = rand_p[i]
        # for G matrix
        self.p_mat = np.full(self.n_snp, -7, dtype=np.float64)  # W = M - 2P
        self.gmat_denom = 0.0  # denominator of G matrix
    # __init__

    def p_matrix(self):
        '''
        P matrix for calculating G matrix
        '''
        for i in range(self.n_snp):
            freq_1 = np.sum(self.par_snp[:, i])/self.n_founder/2.0
            self.p_mat[i] = freq_1
            self.gmat_denom += 2.0*freq_1*(1.0-freq_1)
    # p_matrix

    def msprime_standard(self):
        '''
        To set SNPs in a founder population

        Simulate SNPs under standard Wright-Fisher model.
        call popgen.py
        '''
        self.par_snp = pg.generate_snp(n_sample=self.n_founder,
                                       n_snp=self.n_snp)
    # msprime_standard

    def genotype_matrix(self, pro_inf):
        '''
        Get SNP matrix of progenies

        Args:
            pro_inf (PopulationInfo class): progeny population info
        '''
        pg.progeny_snp(self.par_snp, self.pro_snp, self.snp_dict, pro_inf)
        pg.progeny_genotype(self.pro_snp, self.gen_mat)
    # genotype_matrix
# class SNPInfo


class PhenotypeInfo:
    '''
    class for the information of phenotype

    Phenotypic values,  true/esimated breeding values,
    variance components are stored

    Args:
        mean_phenotype (float): mean phenotype
        residual_var (float): Residial variance

    Attributes:
        mean_pv (float): mean phenotype
        v_e (float): residual variance
        pheno_v (ndarray): phenotypic values
        true_bv (ndarray): true breeding value
        hat_bv (ndarray): estimated breeding value
        hat_beta (ndarray): estimated fixed effects
        hat_vg (float): estimated additvie genetic variance
        hat_ve (float): estimated residual variance
    '''
    def __init__(self, mean_phenotype, residual_var):
        '''
        constructor
        '''
        self.mean_pv = mean_phenotype
        self.v_e = residual_var
        self.pheno_v = []
        self.true_bv = []
        self.hat_bv = []
        self.hat_beta = []
        self.hat_vg = -7.0
        self.hat_ve = -7.0
        # use it to normalize mean true breeding value as 0
        self._first_gen = True
    # __init__

    def calculate_phenotype(self, pro_inf, snp_inf):
        '''
        Calculate phenotype

        Args:
            pro_inf (PopulationInfo class): progeny population
            snp_inf (SNPInfo class): SNP information

        Note:
            The first half of self.pheno_v contains phenotypes of females, and
            seconf half of them are phenotypes of males
        '''
        n_progeny = pro_inf.n_popf + pro_inf.n_popm
        rand_norm = np.random.normal(loc=0.0, scale=np.sqrt(self.v_e),
                                     size=n_progeny)
        self.true_bv = snp_inf.gen_mat @ snp_inf.effect_size.T
        # normalize as average pheno_v is almost equal to mean_pv
        # in the first generation  (F1)
        if self._first_gen:
            bv_mean = np.mean(self.true_bv)
            self.mean_pv -= bv_mean
            self._first_gen = False
        self.pheno_v = self.true_bv + rand_norm + self.mean_pv
    # calclaate_phenotype
# class PhenotypeInfo


def check_tuple(obj, tag, l_tuple):
    '''
    Check an argument

    To check if an argument is tuple and
    if the length of the tuple is correct

    Args:
        obj (unknown): a focal variable
        tag (str): name of the variable
        l_tuple (int): correct tuple length
    '''
    if not isinstance(obj, tuple):
        sys.exit(f'{tag} should be tuple')
    if len(obj) != l_tuple:
        if tag in ('founder_size',  'progeny_size'):
            e_mess = '(No. males, No. females)'
        if tag == 'chrom':
            e_mess = '(chrom no., chrom len (bp), male cM/Mb, female cM/Mb)'
        sys.exit(f'Length of {tag} is incorrect\n{e_mess}')
# check_tuple


class AquaBreeding:
    '''
    class for aquaculture breeding

    Aquaculture breeding is simulated

    Args:
        founder_size (tuple): No. females, no. males (int)
                              in a founder population
        progeny_size (tuple): No. females, no. males (int)
                              in a progeny population
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
                       (int, int, float, float)
        n_snp (int): No. SNPs with effect size
        effect_var (float): Variance of effect size
        mean_phenotype (float): Mean phenotype in F1
        residual_var (float): Redisual variance
        gblup (int): No. SNPs for GBLUP with no effect on phenotype,
                     default None

    Attributes:
        chrom (tuple): chrom num, chrom len, female cM/Mb, male cM/Mb
        par_inf (PopulationInfo class): parental population info
        pro_inf (PopulationInfo class): progeny population info
        snp_inf (SNPInfo class): SNP information with effect size
        phe_inf (PhenotypeInfo class): phenotype, breeding values,
                                       variance components
        cross_inf (ndarray): [[female index, male index, 0, 0],...]
        gblup_inf (SNPInfo class): SNP information for GBLUP
    '''
    def __init__(self, founder_size, progeny_size, chrom, n_snp, effect_var,
                 mean_phenotype, residual_var, gblup=None):
        '''
        constructor
        '''
        # check argument
        check_tuple(founder_size, 'founder_size', 2)
        check_tuple(progeny_size, 'progeny_size', 2)
        check_tuple(chrom, 'chrom', 4)
        # chromosome info
        self.chrom = chrom
        # parental population
        self.par_inf = PopulationInfo(founder_size, self.chrom)
        self.par_inf.init_founder()
        # progeny population
        self.pro_inf = PopulationInfo(progeny_size, self.chrom)
        self.pro_inf.init_progeny()
        # SNP info
        self.snp_inf = SNPInfo(n_snp, effect_var, founder_size, progeny_size,
                               self.chrom)
        # Phenotype info
        self.phe_inf = PhenotypeInfo(mean_phenotype, residual_var)
        # mating info
        self.cross_inf = []
        # SNP for GBLUP
        if gblup is not None:
            self.gblup_inf = SNPInfo(gblup, -1, founder_size, progeny_size,
                                     self.chrom)
        else:
            self.gblup_inf = None
    # __init__

    def print_pop(self, target='founder', n_show=10):
        '''
        Display-print parental/progeny population info

        Args:
            target (str): display 'founder' or 'progeny' population
            n_show (int): No. individuals displayed (female/male each)
        '''
        if target == 'founder':
            self.par_inf.print_popinfo(n_show)
        elif target == 'progeny':
            self.pro_inf.print_popinfo(n_show)
        else:
            sys.exit('target should be founder or progeny')
    # print_pop

    def snp_standard(self):
        '''
        Set SNPs in a founder popolation

        Simulate SNPs using msprime, following standard Wright-Fisher
        model.  Simulate only SNPs with effect size.

        Note:
            Simulate SNPs with demography is under development
        '''
        # SNPs with fixed effect
        self.snp_inf.msprime_standard()
        # SNPs for GBLUP
        if self.gblup_inf is not None:
            self.gblup_inf.msprime_standard()
            self.gblup_inf.p_matrix()
    # snp_standard

    def mating_design(self, cross_design=None, custom_design=None):
        '''
        Define mating design

        Define how the individuals in a founder/parental population
        are crossed.  Call cross_design to use some cross design that is
        already implemented.  Call custom_design to use an user's own mating
        design. custom_design should be numpy.ndarray, where 1st column is
        female parent index, 2nd column is male parent male index, 3rd and 4th
        columns are zero.

        Args:
            cross_design (str): cross design of factorial cross
                                that allows '1x(int)', or 'full' factorial
            custom_design (np.ndarray): female index, male index. no. female
                                        progenies, no. male progenies

        Note:
            * cross_design works when nos. females and males are the same.
            * Evaluate cuastom_design first.  If None, evaluate cross_design.
        '''
        if custom_design is not None:
            self.cross_inf = custom_design
        elif cross_design is not None:
            if self.par_inf.n_popm != self.par_inf.n_popf:
                sys.exit('Nos. males and females shoule be the same')
            self.cross_inf = mt.set_mating_design(self.par_inf, cross_design)
        else:
            sys.exit('Either cross_design or custom_design should be set')
    # mating_design

    def mating(self):
        '''
        Mate founder/parental individuals to produce progenies
        '''
        mt.produce_progeny(self.cross_inf, self.par_inf, self.pro_inf)
    # mating

    def breeding_value(self, blup='ABLUP'):
        '''
        Calculate phenotype and breeding value

        call blup.py

        Args:
            blup (str): If 'ABLUP', numerator relationship matrix is used
                        to estimate breeding values.  If 'GBLUP', genomic
                        relationship matrix is used.  If 'no', breeding
                        values are not estimated.
        '''
        # genotype matrix
        self.snp_inf.genotype_matrix(self.pro_inf)
        # phenotype
        self.phe_inf.calculate_phenotype(self.pro_inf, self.snp_inf)
        # breeding value
        if blup == 'ABLUP':
            bl.ablup(self.phe_inf, self.par_inf, self.pro_inf)
        elif blup == 'GBLUP':
            if self.gblup_inf is None:
                sys.exit('GBLUP option is off.')
            # genotype matrix for GBLUP
            self.gblup_inf.genotype_matrix(self.pro_inf)
            bl.gblup(self.phe_inf, self.gblup_inf)
        elif blup == 'no':
            pass
        else:
            sys.exit('blup should be \'ABLUP\', \'GBLUP\' or \'no\'.')
    # breeding_value

    def selection(self, target='bv', method='mass', top_prop=1.0):
        '''
        Select parents of next generation

        Args:
            target (str): Selection based on breeding value ('bv'),
                          phenotypic value ('phenotype'), or random
                          ('random')
            method (str): How to select from progenies such as mass
                          selection ('mass'), or within-family selection
                          ('within-family')
            top_prop (float): Individuals with top X% breeding values are
                              selected
        '''
        se.select_parent(self.par_inf, self.pro_inf, self.phe_inf,
                         self.cross_inf, target, method, top_prop)
    # selection

    def get_ibd(self):
        '''
        Output inbreeding coefficient in progeny population

        Returns:
            list: inbreeding coefficient
        '''
        ibd_res = []
        # female
        for i in range(self.pro_inf.n_popf):
            tmp_v2 = self.pro_inf.pop_f[i].get_ibd2()
            ibd_res.append(tmp_v2)
        # male
        for i in range(self.pro_inf.n_popm):
            tmp_v2 = self.pro_inf.pop_m[i].get_ibd2()
            ibd_res.append(tmp_v2)
        return ibd_res
    # get_ibd

    def get_phenotype(self):
        '''
        Output phenotype

        Returns:
            list: phenotypics values
        '''
        return list(self.phe_inf.pheno_v)
    # get_phenotype

    def get_true_bv(self):
        '''
        Output true breeding value

        Returns:
            list: true breeding value
        '''
        return list(self.phe_inf.true_bv)
    # get_true_bv

    def get_ebv(self):
        '''
        Output estimated breeding value

        Returns:
            list: estimated breeding value
        '''
        return list(self.phe_inf.hat_bv)
    # get_ebv

    def variance_component(self):
        '''
        Outut additive and residual variance

        Returns:
            tuple: additive and residual variance
        '''
        return (self.phe_inf.hat_vg, self.phe_inf.hat_ve)
    # variance_component
# class AquaBreeding


def main():
    '''
    main
    '''
    print('aquabreeding.py')
    print('Module for simulating aquacuture breeding')
# main


if __name__ == '__main__':
    main()
