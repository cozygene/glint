import replace_missing_values
import logging
from numpy import loadtxt
from tests.test_tools import tools

class MissingValuesTester():
    TEST_042  = "tests/missing_values/files/test_0.42"
    TEST_05 = "tests/missing_values/files/test_0.5"
    TEST_042_05  = "tests/missing_values/files/test_0.42_0.5"
    TEST_STR_IND  = "tests/missing_values/files/test_str_ind"
    TEST_FLOAT_IND  = "tests/missing_values/files/test_float_ind"

    def __init__(self):
        logging.info("Testing Started on MissingValuesTester")
        self.test_replacement()
        logging.info("Testing Finished on MissingValuesTester")

 
    def test_replacement(self):
        logging.info("Start test test_replacement")
        # case 1
        out = replace_missing_values.replace_missing(self.TEST_STR_IND, "NA", 0.42,0.42," ")
        data = loadtxt(self.TEST_042)

        assert data.shape == out.shape
        for i in range(data.shape[0]):
            assert tools.correlation(data[i,:], out[i,:], 1e-2) #-2 because we fist generated data with float64 then we changed it to use float32

        # case 2
        out = replace_missing_values.replace_missing(self.TEST_FLOAT_IND, "9.0", 0.5,0.5," ")
        data = loadtxt(self.TEST_05)

        assert data.shape == out.shape
        for i in range(data.shape[0]):
            assert tools.correlation(data[i,:], out[i,:], 1e-12)

        #case 3
        data = loadtxt(self.TEST_042_05)
        out = replace_missing_values.replace_missing(self.TEST_FLOAT_IND, 9, 0.42,0.5," ")

        assert data.shape == out.shape
        for i in range(data.shape[0]):
            assert tools.correlation(data[i,:], out[i,:], 1e-12)

        logging.info("PASS")

