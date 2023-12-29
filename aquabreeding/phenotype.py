'''
A module for phenotype
'''

import numpy as np
from aquabreeding import blup as bl


def snp_index(n_snp, gblup):
    '''
    Randomly select neutral SNPs and SNPs with effects

    Args:
        n_snp (int): No. causal SNPs
        gblup (int): No. neutral SNPs

    Returns:
        - ndarray: Index of neutral SNPs
        - ndarray: Index of SNPs with effects
    '''
    if gblup is None:
        return None, np.arange(n_snp, dtype=np.int32)
    tmp_index = np.arange(n_snp + gblup, dtype=np.int32)
    np.random.shuffle(tmp_index)
    index_neu = tmp_index[:gblup]
    index_eff = tmp_index[gblup:]
    index_neu.sort()
    index_eff.sort()
    return index_neu, index_eff
# snp_index


def calculate_p_1p(gen_mat, n_hap):
    '''
    Calculate E[p(1 - p)] from a genotype matrix

    Args:
        gen_mat (numpy.ndarray): Genotype matrix (no. samples x np. loci)
        n_hap (int): No. haplotypes

    Returns:
        float: E[p(1 - p)]
    '''
    def cal_freq(arr1d, n_hap):
        '''
        Calculate p(1 - p)

        Args:
            arr1d (numpy.ndarray): List of genotypes
            n_hap (int): No. haplotypes

        Return:
            float: p(1 - p)
        '''
        freq_1 = np.sum(arr1d) / n_hap
        return freq_1 * (1.0 - freq_1)
    return np.mean(np.apply_along_axis(cal_freq, 0, gen_mat, n_hap))
# calculate_p_1p


def effect_size_variance(v_g, gen_eff):
    '''
    Calculate the variance of effect size

    Args:
        v_g (float): Additive genetic variance
        gen_eff (numpy.ndarray): Genotype matrix of the causal SNPs

    Returns:
        float: Variance of effect size
    '''
    n_gen, c_gen = np.shape(gen_eff)
    # E[p(1-p)]
    ep_1p = calculate_p_1p(gen_eff, 2*n_gen)
    # Variance of effect size
    return v_g/c_gen/ep_1p/2.0
# effect_size_variance


def variance_effect(var_p, h2_p, pro_gen):
    '''
    Calculate the variance of effect size

    Args:
        var_p (float): Variance of phenotype
        h2_p (float): Heritability
        pro_gen (numpy.ndarray): Genotype matrix of causal SNPs

    Returns:
        - float: Additive genetic variance
        - float: Residual variance
        - float: Variance of effect size
    '''
    # Additive/residual variance
    v_g = h2_p * var_p
    v_e = (1.0 - h2_p) * var_p
    # Variance of effect size
    v_s = effect_size_variance(v_g, pro_gen)
    return v_g, v_e, v_s
# variance_effect


class PhenotypeInfo:
    '''
    Class for phenotype information

    Phenotypic values,  true/estimated breeding values,
    variance components are stored

    Args:
        mean_phenotype (float): Mean phenotype
        var_phenotype (float): Variance of phenotype
        heritatility (float): Heritability

    Attributes:
        mean_pv (float): Mean phenotype
        var_p (float): Variance of phenotype
        h2_p (float): Heritability
        v_g (float): True additive genetic variance
        v_e (float): True residual variance
        v_s (float): Variance of effect_size
        effect_size (numpy.ndarray): Effect size
        pheno_v (numpy.ndarray): Phenotypic values
        true_bv (numpy.ndarray): True breeding value
        hat_bv (numpy.ndarray): Estimated breeding value
        hat_beta (numpy.ndarray): Estimated fixed effects
        hat_vg (float): Estimated additvie genetic variance
        hat_ve (float): Estimated residual variance
        index_neu (numpy.ndarray): Index of neutral SNPs
        index_eff (numpy.ndarray): Index of causal SNPs
        _first_gen (bool): Check if the fist generation or not
    '''
    def __init__(self, mean_phenotype, var_phenotype, heritability):
        '''
        Constructor
        '''
        self.mean_pv = mean_phenotype
        self.var_p = var_phenotype
        self.h2_p = heritability
        self.v_g = None
        self.v_e = None
        self.v_s = None
        self.effect_size = None
        self.pheno_v = None
        self.true_bv = None
        self.hat_bv = None
        self.hat_beta = None
        self.hat_vg = None
        self.hat_ve = None
        # Index of netural and causal SNPs
        self.index_neu = None
        self.index_eff = None
        # to normalize true breeding value as zero
        self._first_gen = True
    # __init__

    def calculate_bv(self, target, par_inf, pro_inf, founder_size, n_snp,
                     gblup):
        '''
        Calculate phenotype and breeding value

        Args:
            target (str): 'BLUP', 'GBLUP', or 'no'
            par_inf (PopulationInfo): Founder population
            pro_inf (PopulationInfo): Progeny population
            founder_size (tuple): Nos. female and male in the founder
            n_snp (int): No. causal SNPs
            gblup (int): No. neutral SNPs

        Note:
            The first half of self.pheno_v contains phenotypes of females, and
            seconf half of them are phenotypes of males
        '''
        # Set index of neutral and causal SNPs
        if self.index_eff is None:
            self.index_neu, self.index_eff = snp_index(n_snp, gblup)
        # Genotypte matrix of causal SNPs
        gen_eff = pro_inf.gen_mat[:, self.index_eff]
        # Set effect size
        if self.effect_size is None:
            self.v_g, self.v_e, self.v_s = variance_effect(self.var_p,
                                                           self.h2_p,
                                                           gen_eff)
            self.effect_size = np.random.normal(loc=0.0,
                                                scale=np.sqrt(self.v_s),
                                                size=n_snp)
        n_progeny = pro_inf.n_f + pro_inf.n_m
        # true bv
        self.true_bv = gen_eff @ self.effect_size.T
        if self._first_gen:
            self.mean_pv -= np.mean(self.true_bv)
            self._first_gen = False
        # Phenotype
        rand_norm = np.random.normal(loc=0.0, scale=np.sqrt(self.v_e),
                                     size=n_progeny)
        self.pheno_v = self.mean_pv + self.true_bv + rand_norm
        # Numerator relationship matrix
        bl.nrm_cpp(par_inf, pro_inf, founder_size)
        # BLUP
        if target == 'BLUP':
            # Breeding value estimation
            bl.bv_estimation(self, pro_inf.a_mat)
        # GBLUP
        elif target == 'GBLUP':
            # Genomic relationship matrix
            bl.convert_gmatrix(pro_inf, pro_inf.gen_mat[:, self.index_neu])
            # Breeding value estimation
            bl.bv_estimation(self, pro_inf.g_mat)
    # calclaate_phenotype
# PhenotypeInfo


def main():
    '''
    main
    '''
    print('A module for phenotype')
# main


if __name__ == '__main__':
    main()
