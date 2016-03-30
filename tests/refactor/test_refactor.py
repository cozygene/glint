from modules import refactor,methylation_data
from numpy import loadtxt, array_equal
from utils import LinearRegression
import logging
from tests import tools
class FeatureSelectionTester():
    FAKE_DATA  = "tests/refactor/files/feature_selection/data"
    FAKE_CONTROL = "tests/refactor/files/feature_selection/control"
    def __init__(self):
        self.fs_meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA)
        self.test_controls_fs()
        self.test_phenotype_fs()

    def test_controls_fs(self):
        """
        test that all samples that are '1' in the FAKE_CONTROL file are deleted from FAKE_DATA file

        """
        logging.info("Testing controls feature selection...")
        refactor_meth_data = self.fs_meth_data.copy()

        controls = loadtxt(self.FAKE_CONTROL, dtype = str)
        controls_indices = controls[:,1]
        module  = refactor.Refactor(methylation_data = refactor_meth_data, 
                      k = 2, 
                      feature_selection = "controls",
                      phenofile = self.FAKE_CONTROL,
                      t=5)
        index =0
        res = module.feature_selection_handler()

        # controls fs suppose to change methylation data
        assert array_equal(module.meth_data.data, res)

        # res is the data matrix found in FAKE_DATA without the columns (samples) at index i when i is index of a 0 in FAKE_CONTROL file
        for i in range(len(controls_indices)):
            if controls_indices[i]=='0':
                assert array_equal(self.fs_meth_data.data[:,i] ,res[:,index])
                index += 1

        # controls fs result will contain  number of columns (samples) as the number of 0's in the FAKE_CONTROL file
        assert index == len(res[0])

        
        logging.info("PASS")

    def test_phenotype_fs(self):
        logging.info("Testing phenotype feature selection...")
        refactor_meth_data = self.fs_meth_data.copy()

        module  = refactor.Refactor(methylation_data = refactor_meth_data, 
                      k = 2, 
                      feature_selection = "phenotype",
                      phenofile = self.FAKE_CONTROL,
                      t=5)

        phenotype = module._validate_phenotype(self.FAKE_CONTROL, "phenotype")
        # validate phenotype feature selection output (res_data) is correlated to our linear regression for (site, phenotype)
        res_data = module.feature_selection_handler()
        for i,site in enumerate(self.fs_meth_data.data):
            lin_reg = LinearRegression(site, phenotype)
            assert len(lin_reg.residuals) == len(site)
            # validate our residuals are corelated to res_data
            assert tools.correlation(lin_reg.residuals, res_data[i])

        logging.info("PASS")

class RefactorTester():
    DEMO_SMALL_DATA = "tests/refactor/files/demofiles/datafile2"
    DEMO_COVAR = "tests/refactor/files/demofiles/covariates"
    DEMO_PHENO = "tests/refactor/files/demofiles/phenotype"
    DEMO_CELLPRO = "tests/refactor/files/demofiles/cellproportions"

    def __init__(self):
        self.meth_data = methylation_data.MethylationData(datafile = self.DEMO_SMALL_DATA)
        self.test_remove_covariates()
        self.test_low_rank_approx_distances()

    def test_remove_covariates(self):
        logging.info("Testing removing covariates...")
        covar_meth_data = self.meth_data.copy()

        module  = refactor.Refactor(methylation_data = covar_meth_data, 
                      k = 2, 
                      covar = self.DEMO_COVAR,
                      t=500)

        coavr = module._validate_covar(self.DEMO_COVAR,)
        # remove from refactor_meth_data
        module._remove_covariates()

        # remove "manually"
        for i,site in enumerate(self.meth_data.data):
            lin_reg = LinearRegression(site, coavr)
            assert len(lin_reg.residuals) == len(site)

            # validate our residuals are corelated to res_data
            assert tools.correlation(lin_reg.residuals, covar_meth_data.data[i])

        logging.info("PASS")


    def test_low_rank_approx_distances(self):
        """
        tests that number of distances is as the number of sites (distance for every site)
        """
        logging.info("Testing low rank approx distances...")
        dis_meth_data = self.meth_data.copy()

        module  = refactor.Refactor(methylation_data = dis_meth_data, 
                                    k = 5)

        distances = module._calc_low_rank_approx_distances()
        assert distances.size == dis_meth_data.sites_size, "there must be distances as the number of sites"
        logging.info("PASS")

    

    
        