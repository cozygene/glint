from tests.test_tools import tools
from numpy import loadtxt, column_stack
from utils import LinearRegression, LogisticRegression
from modules import methylation_data
import logging

class LinearRegressionTester():
    LIN_REG_X= "tests/utils/files/lin_reg_test_x.txt"
    LIN_REG_Y = "tests/utils/files/lin_reg_test_y.txt"
    LIN_REG_RESIDUALS = "tests/utils/files/lin_reg_test_res.txt"
    LIN_REG_FIT_MODEL = "tests/utils/files/fit_model_results.txt" #In the first row you have the coefficient, f-statistic and p-value of the test on the first column in lin_reg_test_x.txt, and in the second row you have the same details for the second column in lin_reg_test_x.txt
    
    LIN_REG_DATA = "tests/refactor/files/demofiles/datafile2_no_bad_probes"
    LIN_REG_PHENO = "tests/refactor/files/demofiles/phenotype"
    LIN_REG_COVAR = "tests/refactor/files/demofiles/covariates"

    def __init__(self):
        logging.info("Testing Started on LinearRegressionTester")
        self._test_fit_model()
        self._test_regress_out()
        logging.info("Testing Finished on LinearRegressionTester")

    def _test_regress_out(self):
        """
        check linear regression (y-x)
        """
        logging.info("Testing linear regression: regress_out")
        y = loadtxt(self.LIN_REG_Y)
        x = loadtxt(self.LIN_REG_X)
        orig_residuals = loadtxt(self.LIN_REG_RESIDUALS)

        if y.ndim != 1:
            raise("TEST WASNT IMPLEMENTED")
            return

        # test 1 dim
        residuals = LinearRegression.regress_out(y, x)
        y_2dim =  y.reshape(-1, 1)
        residuals2 = LinearRegression.regress_out(y_2dim, x)
        assert tools.correlation(orig_residuals, residuals)
        assert tools.correlation(residuals2, residuals)
        assert len(residuals) == len(x)

        # test 2 dim
        y2 = column_stack((y,y))
        residuals = LinearRegression.regress_out(y2, x)
        for i in range(len(y2[0])):
            assert tools.correlation(orig_residuals, residuals[:,i])
            assert len(residuals[:,i]) == len(x)

        logging.info("PASS")

    def _test_fit_model(self):
        logging.info("Testing linear regression: fit_model")
        meth_data = methylation_data.MethylationDataLoader(datafile = self.LIN_REG_DATA, covarfiles = [self.LIN_REG_COVAR], phenofile = [self.LIN_REG_PHENO])
        results = loadtxt(self.LIN_REG_FIT_MODEL)

        # test 1 dim
        coefs, tstats, pvals = LinearRegression.fit_model(meth_data.phenotype, meth_data.data[0,:], covars = meth_data.covar)
        coefs_inter = coefs[0]
        coefs_site = coefs[-1]
        coefs_covar1 = coefs[1]
        coefs_covar2 = coefs[2]
        tstats = tstats[-1]
        pvals = pvals[-1]

        assert abs(coefs_inter - results[0]) < 1e-3
        assert abs(coefs_site - results[1]) < 1e-3
        assert abs(coefs_covar1 - results[2]) < 1e-3
        assert abs(coefs_covar2 - results[3]) < 1e-3
        assert abs(tstats - results[4]) < 1e-2
        assert abs(pvals - results[5]) < 1e-3
        # Note - there is no option to test 2 dim 
        logging.info("PASS")

class LogisticRegressionTester():
    LIN_REG_DATA = "tests/refactor/files/demofiles/datafile2_no_bad_probes"
    LIN_REG_PHENO = "tests/utils/files/pheno_binary" #binary pheno
    LIN_REG_COVAR = "tests/refactor/files/demofiles/covariates"

    LOG_REG_FIT_MODEL = "tests/utils/files/log_reg_results.txt"

    def __init__(self):
        logging.info("Testing Started on LogisticRegressionTester")
        self.meth_data = methylation_data.MethylationDataLoader(datafile = self.LIN_REG_DATA, covarfiles = [self.LIN_REG_COVAR], phenofile = [self.LIN_REG_PHENO])
        self._test_fit_model()
        logging.info("Testing Finished on LogisticRegressionTester")

    def _test_fit_model(self):
        logging.info("Testing logistic regression: fit_model")
        meth_data = methylation_data.MethylationDataLoader(datafile = self.LIN_REG_DATA, covarfiles = [self.LIN_REG_COVAR], phenofile = [self.LIN_REG_PHENO])
        results = loadtxt(self.LOG_REG_FIT_MODEL)


        coefs, tstats, pvals = LogisticRegression.fit_model(meth_data.phenotype, meth_data.data[0,:], covars = meth_data.covar)
        coefs_inter = coefs[0] 
        coefs_site = coefs[-1]
        coefs_covar1 = coefs[1]
        coefs_covar2 = coefs[2]
        tstats = tstats[-1] #fstat of site under test
        pvals = pvals[-1]
        assert abs(coefs_inter - results[0]) < 1e-3
        assert abs(coefs_site - results[1]) < 1e-3
        assert abs(coefs_covar1 - results[2]) < 1e-3
        assert abs(coefs_covar2 - results[3]) < 1e-3
        assert abs(tstats - results[4]) < 1e-2
        assert abs(pvals - results[5]) < 1e-3
        # Note - there is no option to test 2 dim 
        logging.info("PASS")

