from modules import refactor,methylation_data
from numpy import loadtxt, array_equal
import logging

class FeatureSelectionTester():
    FAKE_DATA  = "tests/refactor/files/feature_selection/data"
    FAKE_CONTROL = "tests/refactor/files/feature_selection/control"
    def __init__(self):
        self.fs_meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA)
        self.test_controls_fs()
        logging.info("PASS")

    def test_controls_fs(self):
        """
        test that all samples that are '1' in the FAKE_CONTROL file are deleted from FAKE_DATA file

        """
        refactor_meth_data = self.fs_meth_data.copy()

        controls = loadtxt(self.FAKE_CONTROL, dtype = str)
        controls_indices = controls[:,1]
        module  = refactor.Refactor(methylation_data = refactor_meth_data, 
                      k = 2, 
                      feature_selection = "controls",
                      phenofile = self.FAKE_CONTROL,
                      t=7)
        index =0
        res = module.feature_selection_handler()

        # res is the data matrix found in FAKE_DATA without the columns (samples) at index i when i is index of a 0 in FAKE_CONTROL file
        for i in range(len(controls_indices)):
            if controls_indices[i]=='0':
                assert array_equal([float(i + x*10) for x in range(10)] ,res[:,index])
                index+=1

        # controls fs result will contain  number of columns (samples) as the number of 0's in the FAKE_CONTROL file
        assert index == len(res[0])

        # controls fs suppose to change methylation data
        assert array_equal(module.meth_data.data, res)

        