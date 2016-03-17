from modules import refactor,methylation_data
from numpy import loadtxt, array_equal, corrcoef
from utils import LinearRegression
import logging

class FeatureSelectionTester():
    FAKE_DATA  = "tests/refactor/files/feature_selection/data"
    FAKE_CONTROL = "tests/refactor/files/feature_selection/control"
    def __init__(self):
        self.fs_meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA)
        self.test_controls_fs()
        self.test_phenotype_fs()
        # logging.info("PASS")

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

        # res is the data matrix found in FAKE_DATA without the columns (samples) at index i when i is index of a 0 in FAKE_CONTROL file
        for i in range(len(controls_indices)):
            if controls_indices[i]=='0':
                assert array_equal([float(i + x*10) for x in range(6)] ,res[:,index])
                index+=1

        # controls fs result will contain  number of columns (samples) as the number of 0's in the FAKE_CONTROL file
        assert index == len(res[0])

        # controls fs suppose to change methylation data
        assert array_equal(module.meth_data.data, res)
        logging.info("PASS")

    def test_phenotype_fs(self):
        logging.info("Testing phenotype feature selection...")
        refactor_meth_data = self.fs_meth_data.copy()

        phenotype = loadtxt(self.FAKE_CONTROL, dtype = str)[:,1].astype(float)
        print phenotype

        module  = refactor.Refactor(methylation_data = refactor_meth_data, 
                      k = 2, 
                      feature_selection = "phenotype",
                      phenofile = self.FAKE_CONTROL,
                      t=5)

        # validate phenotype feature selection output (res_data) is correlated to our linear regression for (site, phenotype)
        res_data = module.feature_selection_handler()
        for i,site in enumerate(refactor_meth_data.data):
            lin_reg = LinearRegression(site, phenotype)
            assert len(lin_reg.residuals) == len(site)
            cor = corrcoef(lin_reg.residuals, res_data[i])

            # validate our residuals are corelated to res_data
            assert abs( 1- cor[0][1]) < 1e-4 

        logging.info("PASS")

    
        