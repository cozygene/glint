from numpy import loadtxt
from utils import tools, LinearRegression
from modules import methylation_data
import logging
from tests.test_tools import tools as tests_tools


class ToolsTester():
    DATA_FILE = "tests/files/datafile2"
    LOW_RANK_APPROX = "tests/utils/files/datafile2_5_rank_approx"
    EUC_DIST = "tests/utils/files/datafile2_euc_dist_first_half_second_half"
    K = 5

    def __init__(self):
        logging.info("Testing Started ToolsTester")
        self.meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA_FILE)
        self.test_low_rank_approx()
        self.test_euclidean_distance()
        logging.info("Testing Finished ToolsTester")

    def test_low_rank_approx(self):
        logging.info("Testing low rank approximation")
        low_rank_approx = loadtxt(self.LOW_RANK_APPROX)
        lra_met_data = self.meth_data.copy()

        res = tools.low_rank_approximation(lra_met_data.data.transpose(), self.K)
        res = res.transpose()

        for i in range(lra_met_data.sites_size):
            assert tests_tools.correlation(res[i,:], low_rank_approx[i,:])
        
        logging.info("PASS")

    def test_euclidean_distance(self):
        logging.info("Testing euclidean distance...")
        euc_dist = loadtxt(self.EUC_DIST)
        out = tools.euclidean_distance(self.meth_data.data[:500,:].transpose(), self.meth_data.data[500:,:].transpose())

        assert tests_tools.correlation(out, euc_dist)
        logging.info("PASS")


class FDRTester():
    DATA = "tests/refactor/files/demofiles/datafile2_no_bad_probes"
    PHENO = "tests/refactor/files/demofiles/phenotype"
    COVAR = "tests/refactor/files/demofiles/covariates"
    QVALUES_RES = "tests/utils/files/qvalues.txt"

    def __init__(self):
        logging.info("Testing Started on WilcoxonTester")
        self.meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA, covarfiles = [self.COVAR], phenofile = self.PHENO)
        self.test_fdr() 
        logging.info("Testing Finished on WilcoxonTester")

    def test_fdr(self):
        logging.info("Testing test_fdr")
        # get pre-calculated qvalues
        qvals_results = loadtxt(self.QVALUES_RES, delimiter=',')

        y = self.meth_data.phenotype # a binary vector (phenotype)
        x = self.meth_data.data.transpose() # all sites  (find pvalues of all sites)

        # glint calc qvalues
        _, _, pvals = LinearRegression.fit_model(y, x, covars = self.meth_data.covar)
        qvals = tools.FDR(pvals)
        # qvals = qvals[2:]

        #check correltaion
        assert(tests_tools.correlation(qvals_results, qvals))
        logging.info("PASS")


class WilcoxonTester():
    DATA = "tests/refactor/files/demofiles/datafile2_no_bad_probes"
    PHENO = "tests/utils/files/pheno_binary" #binary pheno
    COVAR = "tests/refactor/files/demofiles/covariates"

    # U-statistic, p-value results to comare
    U_STATS_RES = 11020.5
    P_VAL_RES = 0.369

    def __init__(self):
        logging.info("Testing Started on WilcoxonTester")
        self.meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA, covarfiles = [self.COVAR], phenofile = self.PHENO)
        self._test_fit_model()
        logging.info("Testing Finished on WilcoxonTester")

    def _test_fit_model(self): # todo not working
        logging.info("Testing Wilcoxon")
        meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA, covarfiles = [self.COVAR], phenofile = self.PHENO)


        y = meth_data.phenotype # a binary vector (phenotype)
        x = meth_data.data[0,:]# site under test - with 0 just the first site

        ustats, pval = tools.wilcoxon_test(y, x)
        import pdb
        print "ustats", ustats
        print "ustats elior", self.U_STATS_RES
        print "pvals", pval
        print "pval elior", self.P_VAL_RES
        pdb.set_trace() #somethin is wrong here data - all or just one site????
        assert abs(ustats - self.U_STATS_RES) < 1e-3 #todo something is wrong
        assert abs(pval - self.P_VAL_RES) < 1e-3
        logging.info("PASS")