import replace_missing_values
import logging
from numpy import loadtxt

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
        out = replace_missing_values.replace_missing(self.TEST_STR_IND, "NA", 0.42,0.42," ")
        data = loadtxt(self.TEST_042)
        assert (out == data).all()

        out = replace_missing_values.replace_missing(self.TEST_FLOAT_IND, "9.0", 0.5,0.5," ")
        data = loadtxt(self.TEST_05)
        assert (out == data).all()


        data = loadtxt(self.TEST_042_05)
        out = replace_missing_values.replace_missing(self.TEST_FLOAT_IND, 9, 0.42,0.5," ")
        assert (out == data).all()
        logging.info("PASS")

