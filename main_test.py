import validate_dependencies
from configuration import configurelogging
import logging
configurelogging.ConfigureLogging(loglevel = logging.DEBUG)
from tests.refactor import test_refactor
from tests.ewas import test_ewas
from tests.utils import test_lin_reg, test_pca, test_tools
from tests.methylation_data import test_methylation_data
from tests.imputing import test_imputing
from tests.missing_values import replace_missing_values_test
from tests.houseman import test_houseman

test_refactor.RefactorTester()
test_refactor.SenariosTester()
test_lin_reg.LinearRegressionTester()
test_lin_reg.LogisticRegressionTester()
test_methylation_data.DataTester()
test_pca.PCATester()
test_tools.ToolsTester()
test_tools.WilcoxonTester()
test_tools.FDRTester()
# test_kit.PCAKitTester()
test_imputing.ImputingTester()
replace_missing_values_test.MissingValuesTester()
test_ewas.LMMTester()
test_houseman.HousemanTester()