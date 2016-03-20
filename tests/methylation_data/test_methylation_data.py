from modules import methylation_data
from numpy import loadtxt, array_equal
import logging

class DataTester():
    FAKE_DATA  = "tests/methylation_data/files/std/data.txt"
    FAKE_DATA_STDTH =  "tests/methylation_data/files/std/data_after_std_0.25.txt"
    STDTH = 0.25

    def __init__(self):
        self.meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA)
        self.test_stdth()

    def test_stdth(self):
        logging.info("Testing stdth...")
        data_copy = self.meth_data.copy()
        data_copy.remove_lowest_std_sites(self.STDTH)
        data_after_std = methylation_data.MethylationData(datafile = self.FAKE_DATA_STDTH)
        assert array_equal(data_copy.data, data_after_std.data)
        logging.info("PASS")