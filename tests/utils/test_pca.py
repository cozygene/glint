from numpy import loadtxt
from utils import pca
from modules import methylation_data
import logging
from tests.test_tools import tools

class PCATester():
    DATA_FILE = "tests/files/datafile2"
    PCA_P_RES = "tests/utils/files/datafile2_pcs"

    def __init__(self):
        logging.info("Testing Started on PCATester")
        pca_res_p = loadtxt(self.PCA_P_RES)

        meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA_FILE)
        pca_out = pca.PCA(meth_data.data.transpose())

        for i in range(10):
            assert tools.correlation(pca_out.P[:,i], pca_res_p[:,i])

        logging.info("PASS")
        logging.info("Testing Finished on PCATester")