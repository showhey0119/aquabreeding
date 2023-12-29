'''
A module for population, individual, and chromosome class
'''

import numpy as np
from aquabreeding import aquacpp as cpp


def merge_snp(snp_pro, pop_x, stt):
    '''
    Generate SNP matrix for female or male

    Args:
        snp_pro (numpy.ndarray): SNP matrix
        pop_x (list): List of IndividualInfo of female/male
        stt (int): Start index of rows
    '''
    i_row = stt
    for individual in pop_x:
        i_col = 0
        for chrom in individual.chrom_ls:
            n_pos = np.shape(chrom.position)[0]
            cpp.copy_1D_to_2D(chrom.snp_mat, snp_pro, i_row, i_col)
            cpp.copy_1D_to_2D(chrom.snp_pat, snp_pro, i_row+1, i_col)
            i_col += n_pos
        i_row += 2
# merge_snp


class ChromInfo:
    '''
    Class for a chromosome

    Args:
        chrom (tuple): Chrom num, chrom len, female cM/Mb, male cM/Mb

    Attributes:
        chrom (tuple): Chrom num, chrom len, female cM/Mb, male cM/Mb
        position (numpy.ndarray): Positions of SNPs
        snp_mat (numpy.ndarray): SNPs of the maternal chromosome
        snp_pat (numpy.ndarray): SNPs of the paternal chromosome
    '''
    def __init__(self, chrom):
        '''
        Constractor
        '''
        self.chrom = chrom
        # position of SNPs
        self.position = None
        # SNPs in mat/paternal chromosomes
        self.snp_mat = None
        self.snp_pat = None
    # __init__
# ChromInfo


class IndividualInfo:
    '''
    Class for an individual and its genome

    Args:
        chrom (tuple): Chrom num, chrom len, female cM/Mb, male cM/Mb
        sample_id (int): ID of an individual, default = -1
        pat_id (int): ID of father, default = -1
        mat_id (int): ID of mother, defautl = -1

    Attributes:
        sample_id (int): ID of the parents
        chrom (tuple): Chrom num, chrom len, female cM/Mb, male cM/Mb
        chrom_ls (list): List of ChromInfo class
        mat_id (int): ID of mother
        pat_id (int): ID of father
        n_recf (float): Expected num crossing-over event in female
        n_recm (float): Expected num crossing-over event in male
    '''
    def __init__(self, chrom, sample_id=-1, pat_id=-1, mat_id=-1):
        '''
        Constructor
        '''
        self.sample_id = sample_id
        self.chrom = chrom
        # paternal and maternal genomes
        self.chrom_ls = []
        for _ in range(self.chrom[0]):
            self.chrom_ls.append(ChromInfo(self.chrom))
        # ID of parents (-1 means unknown)
        self.pat_id = pat_id
        self.mat_id = mat_id
        # expected num crossing-over event
        self.n_recf = 2.0*self.chrom[2]*self.chrom[1]/(10.0**6)/100.0
        self.n_recm = 2.0*self.chrom[3]*self.chrom[1]/(10.0**6)/100.0
    # __init__
# IndividualInfo


class PopulationInfo:
    '''
    Class for breeding population information

    Args:
        pop_size (tuple): Nos. females and  males in a population
        chrom (tuple): Chrom num, chrom len, female cM/Mb, male cM/Mb
        n_nat (int): No. additional wild individuals

    Attributes:
        chrom (tuple): Chrom num, chrom len, female cM/Mb, male cM/Mb
        n_f (int): No. females in a population
        n_m (int): No. males in a population
        pop_f (list): List of female IndividualInfo class
        pop_m (list): List of male IndividualInfo class
        tmp_id (int): ID of all parents (only used for a founder)
        d_par (dict): Dict of parents of founders (only used for a founder)
        n_nat (int): No. additional wild individual
        pop_n (list): List of additional wild IndividualInfo class
        count_nat (int): Counter for additional wilds
        gen_mat (numpy.ndarray): Genotype matrix
        a_mat (numpy.ndarray): Numerator relationship matrix
        g_mat (numpy.ndarray): Genomic relationship matrix
    '''
    def __init__(self, pop_size, chrom, n_nat=None):
        '''
        Constructur
        '''
        self.chrom = chrom
        # No. individuals
        self.n_f = pop_size[0]
        self.n_m = pop_size[1]
        # List of individual class
        self.pop_f = []
        self.pop_m = []
        # Pedigree info
        self.tmp_id = 0
        self.d_par = {}
        # For additional wild individuals
        self.n_nat = n_nat
        self.pop_n = []
        self.count_nat = 0
        # Genotype matrix
        self.gen_mat = None
        # Numerator/genomic relationship matrix
        self.a_mat = None
        self.g_mat = None
    # __init__

    def init_founder(self):
        '''
        Initiate founder/parental population
        '''
        # female
        for _ in range(self.n_f):
            self.pop_f.append(IndividualInfo(self.chrom, self.tmp_id))
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = -1
            self.d_par[self.tmp_id]['pat'] = -1
            self.tmp_id += 1
        # male
        for _ in range(self.n_m):
            self.pop_m.append(IndividualInfo(self.chrom, self.tmp_id))
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = -1
            self.d_par[self.tmp_id]['pat'] = -1
            self.tmp_id += 1
        # extral natural individuals
        if self.n_nat is not None:
            for _ in range(self.n_nat):
                self.pop_n.append(IndividualInfo(self.chrom, self.tmp_id))
                self.d_par[self.tmp_id] = {}
                self.d_par[self.tmp_id]['mat'] = -1
                self.d_par[self.tmp_id]['pat'] = -1
                self.tmp_id += 1
    # init_founder

    def init_progeny(self):
        '''
        Initiate progeny population
        '''
        # female
        for _ in range(self.n_f):
            self.pop_f.append(IndividualInfo(self.chrom))
        # male
        for _ in range(self.n_m):
            self.pop_m.append(IndividualInfo(self.chrom))
    # init_progeny

    def change_size(self, new_size):
        '''
        Change nos. females and males

        Args:
            new_size (tuple): New nos. females and males
        '''
        # female
        if self.n_f > new_size[0]:
            del self.pop_f[new_size[0]:]
        elif self.n_f < new_size[0]:
            for _ in range(new_size[0] - self.n_f):
                self.pop_f.append(IndividualInfo(self.chrom))
        self.n_f = new_size[0]
        # male
        if self.n_m > new_size[1]:
            del self.pop_m[new_size[1]:]
        elif self.n_m < new_size[1]:
            for _ in range(new_size[1] - self.n_m):
                self.pop_m.append(IndividualInfo(self.chrom))
        self.n_m = new_size[1]
    # change_size

    def genotype_matrix(self, n_snp, gblup):
        '''
        Generate the genotype matrix

        Args:
            n_snp (int): No. causal SNPs
            gblup (int): No. neutral SNPs
        '''
        n_total = self.n_f + self.n_m
        if gblup is None:
            all_snp = n_snp
        else:
            all_snp = n_snp + gblup
        snp_pro = np.full((2*n_total, all_snp), -7, dtype=np.int32)
        # female
        merge_snp(snp_pro, self.pop_f, 0)
        # male
        merge_snp(snp_pro, self.pop_m, 2*self.n_f)
        # Convert SNP to genotype matrix
        self.gen_mat = np.full((n_total, all_snp), -7, dtype=np.int32)
        cpp.snp_to_genotype(self.gen_mat, snp_pro)
    # genotype_matrix

    def new_founder_id(self):
        '''
        Set IDs for new founder
        '''
        # female
        for individual in self.pop_f:
            individual.sample_id = self.tmp_id
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = individual.mat_id
            self.d_par[self.tmp_id]['pat'] = individual.pat_id
            self.tmp_id += 1
        # male
        for individual in self.pop_m:
            individual.sample_id = self.tmp_id
            self.d_par[self.tmp_id] = {}
            self.d_par[self.tmp_id]['mat'] = individual.mat_id
            self.d_par[self.tmp_id]['pat'] = individual.pat_id
            self.tmp_id += 1
    # new_founder_id
# PopulationInfo class


def main():
    '''
    main
    '''
    print('A module for population, individual, and chromosome class')
# main


if __name__ == '__main__':
    main()
