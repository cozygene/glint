from configuration import configurelogging
configurelogging.configureLogging("test")
from tests.refactor import test_refactor
from tests.ewas import test_ewas
from tests.utils import test_lin_reg, test_pca
from tests.methylation_data import test_methylation_data
test_refactor.RefactorTester()
test_refactor.FeatureSelectionTester()
test_lin_reg.LinearRegressionTester()
test_methylation_data.DataTester()
test_pca.PCATester()
