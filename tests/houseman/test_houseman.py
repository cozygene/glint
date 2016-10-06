from modules import houseman, methylation_data
from parsers import houseman_parser
from numpy import loadtxt
import logging
from tests.test_tools import tools

class HousemanTester():
    DATA  = "tests/houseman/files/datafile.txt"
    OUT_COMP  = "tests/houseman/files/datafile_houseman_est_matlab.txt"

    def __init__(self):
        logging.info("Testing Started on HousemanTester")
        self.meth_data = methylation_data.MethylationDataLoader(datafile = self.DATA)
        self.test_components()
        logging.info("Testing Finished on HousemanTester")

    def test_components(self):
        logging.info("Testing houseman components...")

        comp = loadtxt(self.OUT_COMP)

        module  = houseman.Houseman(self.meth_data, open(houseman_parser.HOUSEMAN_DEFAULT_REFERENCE, 'r'))
        
        module.components
        assert module.components.shape == comp.shape
        for i in range(module.components.shape[1]):
            assert tools.correlation(module.components[:,i], comp[:,i], 1e-2)
        logging.info("PASS")