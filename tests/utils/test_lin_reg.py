from tests.test_tools import tools
from numpy import loadtxt, column_stack
from utils import LinearRegression
import logging

class LinearRegressionTester():
    LIN_REG_X= "tests/utils/files/lin_reg_test_x.txt"
    LIN_REG_Y = "tests/utils/files/lin_reg_test_y.txt"
    LIN_REG_RESIDUALS = "tests/utils/files/lin_reg_test_res.txt"
    LIN_REG_FIT_MODEL = "tests/utils/files/fit_model_results.txt" #In the first row you have the coefficient, f-statistic and p-value of the test on the first column in lin_reg_test_x.txt, and in the second row you have the same details for the second column in lin_reg_test_x.txt
    
    def __init__(self):
        self.y = loadtxt(self.LIN_REG_Y)
        self.x = loadtxt(self.LIN_REG_X)
        self._test_fit_model()
        self._test_regress_out()

    def _test_regress_out(self):
        """
        check linear regression (y-x)
        """
        logging.info("Testing linear regression: regress_out")
        orig_residuals = loadtxt(self.LIN_REG_RESIDUALS)

        if self.y.ndim != 1:
            raise("TEST WASNT IMPLEMENTED")
            return

        # test 1 dim
        residuals = LinearRegression.regress_out(self.y, self.x)
        y_2dim =  self.y.reshape(-1, 1)
        residuals2 = LinearRegression.regress_out(y_2dim, self.x)
        assert tools.correlation(orig_residuals, residuals)
        assert tools.correlation(residuals2, residuals)
        assert len(residuals) == len(self.x)

        # test 2 dim
        y2 = column_stack((self.y,self.y))
        residuals = LinearRegression.regress_out(y2, self.x)
        for i in range(len(y2[0])):
            assert tools.correlation(orig_residuals, residuals[:,i])
            assert len(residuals[:,i]) == len(self.x)

        logging.info("PASS")

    def _test_fit_model(self):
        logging.info("Testing linear regression: fit_model")
        results = loadtxt(self.LIN_REG_FIT_MODEL)
        # test 1 dim
        for i in range(self.x.shape[1]):
            coefs, fstats, pvals = LinearRegression.fit_model(self.y, self.x[:,i])
            assert abs(coefs - results[i][0]) < 1e-3
            assert abs(fstats - results[i][1]) < 1e-3
            assert abs(pvals - results[i][2]) < 1e-3
        
        # test 2 dim 
        coefs, fstats, pvals = LinearRegression.fit_model(self.y, self.x)
        for i in range(self.x.shape[1]):
            assert abs(coefs[i] - results[i][0]) < 1e-3
            assert abs(fstats[i] - results[i][1]) < 1e-3
            assert abs(pvals[i] - results[i][2]) < 1e-3
        logging.info("PASS")
