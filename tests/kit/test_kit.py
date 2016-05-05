from modules import methylation_data, kit
import logging
 
class PCAKitTester(object):
    FAKE_DATA  = "tests/kit/files/data.txt"
    def __init__(self):
        self.meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA)
        self.pca_kit = kit.PCAKit(self.meth_data.copy())

        self.test_remove_outliers()

    def test_remove_outliers(self):
        logging.info("Test remove outliers")
        self.pca_kit.draw_pca_scatter(3)
        # raw_input()
        self.pca_kit.exclude_maxpcstds([[3,1]])
        self.pca_kit.draw_pca_scatter(5)
        logging.info("PASS")