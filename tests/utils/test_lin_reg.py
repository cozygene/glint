
from numpy import loadtxt, corrcoef, empty_like
from utils import tools, pca, LinearRegression, common
import logging

class LinearRegressionTester():
    LIN_REG_X= "tests/utils/files/lin_reg_test_x.txt"
    LIN_REG_Y = "tests/utils/files/lin_reg_test_y.txt"
    LIN_REG_RESIDUALS = "tests/utils/files/lin_reg_test_res.txt"

    def __init__(self):
        """
        check linear regression (y-x)
        """
        y = loadtxt(self.LIN_REG_Y)
        x = loadtxt(self.LIN_REG_X)
        residuals = loadtxt(self.LIN_REG_RESIDUALS)

         
        if len(y.shape) == 1:
            lin_reg = LinearRegression(y, x)
            cor = corrcoef(lin_reg.residuals, residuals)
            assert abs( 1- cor[0][1]) < 1e-4 
            assert len(lin_reg.residuals) == len(x)

        elif len(y.shape) == 2:  
            raise("TEST WASNT IMPLEMENTED")

        logging.info("PASS")