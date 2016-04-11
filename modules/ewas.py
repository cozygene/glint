from utils import LinearRegression #, plot
from numpy import column_stack
from module import Module
from utils import common
import logging

class EWAS(Module):
    AVALIABLE_TESTS = ['linear_regression', 'logistic_regression']
    TEST_FUNC_NAME_FORMAT = "_{test_name}_test"   # feature selections function name format


    def __init__(self, methylation_data, test):
        self.meth_data = methylation_data
        self.test_handlers = self._get_test_handler(test)

    def run(self):
        logging.info('starting EWAS...');
        [test_handler() for test_handler in self.test_handlers]
        logging.info('EWAS is Done!')

    def _get_test_handler(self, tests_list):
        # check that the tests in test_list are all optional tests (found in AVALIABLE_TESTS)
        if set(set(tests_list).difference(set(self.AVALIABLE_TESTS))) == 0:
            common.terminate('tests %s are not available' % str(set(tests_list).difference(set(self.AVALIABLE_TESTS))))

        return [getattr(self, self.TEST_FUNC_NAME_FORMAT.format(test_name=test)) for test in tests_list]

    def _logistic_regression_test(self):
        pass
        
    def _linear_regression_test(self):
        """
        linear regression test
        """
        print self.meth_data.phenotype.shape
        if self.meth_data.covar:
            test = column_stack((self.meth_data.phenotype, self.meth_data.covar))
        else:
            test = self.meth_data.phenotype.copy()
        print test.shape
        lin_reg = LinearRegression(self.meth_data.data, test)
        
        # result is list of tuples (p_value, t_statistix, coef) for each site
        # or maybe (cpgname, p_value, t_statistix, coef)

    def output_result():
        #save to file a table of the format
        # cpgname | t_statistix | coef | p_value
        # sorted by p_value (will not keep the order of cpgname) when lowest p_value will be at the top of the table
        pass

    def save_qqplot():
        pass
        # plot.draw_qqplot(y= p_values, title='bla bla TODO change this', xtitle='-log10(expected)', ytitle='-log10(observed)')