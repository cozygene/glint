
from numpy import loadtxt, corrcoef
from utils import tools, pca, LinearRegression, common
import logging

class LinearRegressionTester():
    LIN_REG_DATA = "tests/ewas/files/lin_reg_test_data.txt"
    LIN_REG_PHENO = "tests/ewas/files/lin_reg_test_pheno.txt"
    LIN_REG_RESIDUALS = "tests/ewas/files/lin_reg_test_res.txt"

    def __init__(self):
        data = loadtxt(self.LIN_REG_DATA)
        pheno = loadtxt(self.LIN_REG_PHENO)
        residuals = loadtxt(self.LIN_REG_RESIDUALS)

        lin_reg = LinearRegression(pheno, data)

        assert len(lin_reg.residuals) == len(pheno)

        cor = corrcoef(lin_reg.residuals, residuals)

        # validate our residuals are corelated to LIN_REG_RESIDUALS
        assert abs( 1- cor[0][1]) < 1e-4 

        logging.info("PASS")