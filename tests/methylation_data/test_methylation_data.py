from modules import methylation_data
from numpy import loadtxt, array_equal
import logging
from tests import tools

class DataTester():
    FAKE_DATA  = "tests/methylation_data/files/data.txt"
    FAKE_PHENO = "tests/methylation_data/files/pheno.txt"
    FAKE_COVAR = "tests/methylation_data/files/covar.txt"
    FAKE_DATA_STDTH = "tests/methylation_data/files/data_after_std_0.25.txt"
    FAKE_DATA_INC = "tests/methylation_data/files/data_inc_1_2_3.txt"
    FAKE_DATA_EXC = "tests/methylation_data/files/data_ex_1_2_3.txt"
    FAKE_DATA_KEEP = "tests/methylation_data/files/data_keep_1_2_3.txt"
    FAKE_DATA_REMOVE = "tests/methylation_data/files/data_remove_1_2_3.txt"
    FAKE_PHENO_KEEP = "tests/methylation_data/files/pheno_keep_1_2_3.txt"
    FAKE_PHENO_REMOVE = "tests/methylation_data/files/pheno_remove_1_2_3.txt"
    FAKE_COVAR_KEEP = "tests/methylation_data/files/covar_keep_1_2_3.txt"
    FAKE_COVAR_REMOVE = "tests/methylation_data/files/covar_remove_1_2_3.txt"
    FAKE_DATA_MAX_MEANS = "tests/methylation_data/files/data_excluded_mean_above_0.5.txt"
    FAKE_DATA_MIN_MEANS = "tests/methylation_data/files/data_excluded_mean_below_0.5.txt"
    FAKE_DATA_MEANS= "tests/methylation_data/files/data_means.txt"
    MIN_MEAN_TH = 0.5
    MAX_MEAN_TH = 0.5
    STDTH = 0.25
    INC_EXC = ['cg1', 'cg2', 'cg3']
    KEEP_REMOVE_INDICES = ['sample1', 'sample2', 'sample3']

    def __init__(self):
        self.meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA, covarfiles = [self.FAKE_COVAR], phenofile = self.FAKE_PHENO)
        self.test_remove_lowest_std_sites()
        self.test_get_mean_per_site()
        self.test_include()
        self.test_exclude()
        self.test_keep()
        self.test_remove()
        self.test_exclude_sites_with_low_mean()
        self.test_exclude_sites_with_high_mean()
        self.test_upload_new_files()


    def test_remove_lowest_std_sites(self):
        logging.info("Testing stdth...")
        data_copy = self.meth_data.copy()
        data_copy.remove_lowest_std_sites(self.STDTH)
        data_after_std = methylation_data.MethylationData(datafile = self.FAKE_DATA_STDTH)
        assert array_equal(data_copy.data, data_after_std.data)
        logging.info("PASS")


    def test_include(self):
        logging.info("Testing include...")
        data_after = methylation_data.MethylationData(datafile = self.FAKE_DATA_INC)
        data = self.meth_data.copy()
        data.include(self.INC_EXC)
        assert array_equal(data_after.data, data.data)
        logging.info("PASS")

    def test_exclude(self):
        logging.info("Testing test_exclude...")
        data_after = methylation_data.MethylationData(datafile = self.FAKE_DATA_EXC)
        data = self.meth_data.copy()
        data.exclude(self.INC_EXC)
        assert array_equal(data_after.data, data.data)
        logging.info("PASS")
        
    def test_keep(self):
        logging.info("Testing keep...")
        data_after = methylation_data.MethylationData(datafile = self.FAKE_DATA_KEEP, covarfiles = [self.FAKE_COVAR_KEEP], phenofile = self.FAKE_PHENO_KEEP)
        data = self.meth_data.copy()
        data.keep(self.KEEP_REMOVE_INDICES)
        assert array_equal(data_after.data, data.data)
        assert array_equal(data_after.phenotype, data.phenotype)
        assert array_equal(data_after.covar, data.covar)
        logging.info("PASS")

    def test_remove(self):
        logging.info("Testing remove...")
        data_after = methylation_data.MethylationData(datafile = self.FAKE_DATA_REMOVE, covarfiles = [self.FAKE_COVAR_REMOVE], phenofile = self.FAKE_PHENO_REMOVE)
        data = self.meth_data.copy()
        data.remove(self.KEEP_REMOVE_INDICES)
        assert array_equal(data_after.data, data.data)
        assert array_equal(data_after.phenotype, data.phenotype)
        assert array_equal(data_after.covar, data.covar)
        logging.info("PASS")

    def test_get_mean_per_site(self):
        logging.info("Testing mean oer site...")
        meth_data_means = self.meth_data.get_mean_per_site()
        means = loadtxt(self.FAKE_DATA_MEANS)
        assert tools.correlation(means, meth_data_means)
        logging.info("PASS")

    def test_exclude_sites_with_low_mean(self):
        logging.info("Testing excluded mean below 0.5...")
        data = self.meth_data.copy()
        data.exclude_sites_with_low_mean(self.MIN_MEAN_TH)
        res = loadtxt(self.FAKE_DATA_MIN_MEANS)
        assert array_equal(res, data.data)
        logging.info("PASS")

    def test_exclude_sites_with_high_mean(self):
        logging.info("Testing excluded mean above 0.5...")
        data = self.meth_data.copy()
        data.exclude_sites_with_high_mean(self.MAX_MEAN_TH)
        res = loadtxt(self.FAKE_DATA_MAX_MEANS)
        assert array_equal(res, data.data)
        logging.info("PASS")

    def test_upload_new_files(self):
        logging.info("Testing upload new covaritates and phenotype files...")
        data = self.meth_data.copy()
        data_upload = methylation_data.MethylationData(datafile = self.FAKE_DATA_REMOVE)
        
        data.remove(self.KEEP_REMOVE_INDICES)

        data_upload.upload_new_covaritates_files([self.FAKE_COVAR_REMOVE])
        data_upload.upload_new_phenotype_file(self.FAKE_PHENO_REMOVE)


        assert array_equal(data.data, data_upload.data)
        assert array_equal(data.phenotype, data_upload.phenotype)
        assert array_equal(data.covar, data_upload.covar)
        logging.info("PASS")
