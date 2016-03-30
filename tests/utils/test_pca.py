from numpy import loadtxt
from utils import pca
from modules import methylation_data
import logging
from tests import tools

class PCATester():
    DATA_FILE = "tests/utils/files/datafile2"
    PCA_P_RES = "tests/utils/files/datafile2_pcs"

    def __init__(self):
        pca_res_p = loadtxt(self.PCA_P_RES)

        meth_data = methylation_data.MethylationData(datafile = self.DATA_FILE)
        pca_out = pca.PCA(meth_data.data.transpose())

        for i in range(10):
            assert tools.correlation(pca_out.P[:,i], pca_res_p[:,i])

        logging.info("PASS")