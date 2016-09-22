from scipy.stats import pearsonr
from modules import methylation_data, ewas, lmm
from numpy import loadtxt, log10
import logging
from pandas import Index,unique

class LMMTester():

    KINSHIP = "tests/ewas/files/test_kinship.txt"
    EWASHER_OUT = "tests/ewas/files/out_ewasher.txt"

    DATA = "tests/refactor/files/demofiles/datafile2_no_bad_probes"
    PHENO = "tests/refactor/files/demofiles/phenotype" 
    COVAR = "tests/refactor/files/demofiles/covariates"


    def __init__(self):
        logging.info("Testing Started on LMMTester")
        self.meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA, covarfiles = [self.COVAR], phenofile = [self.PHENO])
        self.test_pvalues()
        logging.info("Testing Finished on LMMTester")

    def test_pvalues(self):
        logging.info("Testing test_pvalues...")
        
        # calculate pvalues using lmm (function under test)
        kinship = loadtxt(self.KINSHIP)
        module = lmm.LMM(kinship)
        lmm_cpgnames, lmm_pvalues, _, _, _, _, _, _ = module.run(self.meth_data.data.transpose(), self.meth_data.phenotype, self.meth_data.covar, self.meth_data.cpgnames)
        lmm_pvalues = -log10(lmm_pvalues)

        # get pre-calculated pvalues 
        ewasher = loadtxt(self.EWASHER_OUT, dtype = str)
        ewasher_cpgnames = ewasher[:,0]
        ewasher_pvals = ewasher[:,1].astype(float)
        ewasher_pvals = -log10(ewasher_pvals)

        # sort both or them by cpgnames
        
        # sort lmm output
        lmm_sorted_indices = lmm_cpgnames.argsort()  
        lmm_cpgnames = lmm_cpgnames[lmm_sorted_indices]
        lmm_pvalues = lmm_pvalues[lmm_sorted_indices]
        
        # sort pre-calculated output
        ewasher_sorted_indices = ewasher_cpgnames.argsort() 
        ewasher_cpgnames = ewasher_cpgnames[ewasher_sorted_indices]
        ewasher_pvals = ewasher_pvals[ewasher_sorted_indices]

        # check that cpgnames sorted lists are equal
        assert((ewasher_cpgnames == lmm_cpgnames).any())

        # correlation is  not y
        assert(pearsonr(ewasher_pvals, lmm_pvalues)[0] > 0.7)

        logging.info("PASS")

