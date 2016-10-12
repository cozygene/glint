from modules import imputing,methylation_data
from numpy import loadtxt, array_equal
import logging
from tests.test_tools import tools

class ImputingTester():

    SNPS_FILE = "tests/imputing/files/reut22_5snps.snps"
    GENO_FILE = "tests/imputing/files/reut22_10samples_5snps.geno"
    MISSING_SAMPLES_GENO_FILE = "tests/imputing/files/reut22_10samples_5snps-missing_samples.geno"
    MISSING_SNPS_GENO_FILE = "tests/imputing/files/reut22_10samples_5snps-missing_snps.geno"
    MISSING_SAMPLES_AND_SNPS_GENO_FILE = "tests/imputing/files/reut22_10samples_5snps-missing_samples_and_snps.geno"
    IND_FILE = "tests/imputing/files/reut22_10s.ind"
    MANUALLY_IMPUTED_FILE = "tests/imputing/files/reut22_test_out" 
    MISSING_SAMPLES_MANUALLY_IMPUTED_FILE = "tests/imputing/files/reut22_test_out-missing_samples" 
    MISSING_SNPS_MANUALLY_IMPUTED_FILE = "tests/imputing/files/reut22_test_out-missing_snps" 
    MISSING_SAMPLES_AND_SNPS_MANUALLY_IMPUTED_FILE = "tests/imputing/files/reut22_test_out-missing_samples_and_snps"
    SITES_SCORES_FILE = 'parsers/assets/sites_scores_list'
    SITES_SNPS_FILE = 'parsers/assets/site_snps_list'
    SITES_IDS_FILE ='parsers/assets/sites_ids_list'
    SNPS_IDS_FILE = 'parsers/assets/snps_ids_list'
    SITES_SNPS_COEFF_FILE = 'parsers/assets/sites_snps_coeff_list'

    MIN_SCORE = 0.0005
    MIN_MISSING_VALUES = 0.2
    def __init__(self):
        logging.info("Testing Started on ImputingTester")
        self.module  = imputing.Imputation(self.SITES_SCORES_FILE, self.SITES_SNPS_FILE, self.SITES_IDS_FILE, self.SNPS_IDS_FILE, self.SITES_SNPS_COEFF_FILE)
        self.test_impute()
        # self.test_impute_missing_samples() #missing samples handling
        self.test_impute_missing_snps()
        # self.test_impute_missing_samples_and_snps() #missing samples handling
        logging.info("Testing Finished on ImputingTester")

    def test_impute(self):
        logging.info("Testing test_impute...")
        imputed = loadtxt(self.MANUALLY_IMPUTED_FILE)
        self.module.impute(self.MIN_SCORE, self.SNPS_FILE, self.GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        
        assert(imputed.shape[0] == self.module.site_imputation.shape[0])
        assert(imputed.shape[1] == self.module.site_imputation.shape[1])
        for i in range(imputed.shape[1]):
            assert(tools.correlation(imputed[:,i], self.module.site_imputation[:,i]))
        
        for i in range(imputed.shape[0]):
            assert(tools.correlation(imputed[i,:], self.module.site_imputation[i,:]))
        logging.info("PASS")

    # def test_impute_missing_samples(self):  #missing samples handling
    #     logging.info("Testing test_impute_missing_samples...")
    #     imputed = loadtxt(self.MISSING_SAMPLES_MANUALLY_IMPUTED_FILE)
    #     self.module.impute(self.MIN_SCORE, self.SNPS_FILE, self.MISSING_SAMPLES_GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
    #     assert(imputed.shape[0] == self.module.site_imputation.shape[0])
    #     assert(imputed.shape[1] == self.module.site_imputation.shape[1])
    #     for i in range(imputed.shape[1]):
    #         assert(tools.correlation(imputed[:,i], self.module.site_imputation[:,i]))
        
    #     for i in range(imputed.shape[0]):
    #         assert(tools.correlation(imputed[i,:], self.module.site_imputation[i,:]))
    #     logging.info("PASS")


    def test_impute_missing_snps(self):
        logging.info("Testing test_impute_missing_snps...")
        imputed = loadtxt(self.MISSING_SNPS_MANUALLY_IMPUTED_FILE)
        self.module.impute(self.MIN_SCORE, self.SNPS_FILE, self.MISSING_SNPS_GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        
        assert(imputed.shape[0] == self.module.site_imputation.shape[0])
        assert(imputed.shape[1] == self.module.site_imputation.shape[1])

        for i in range(imputed.shape[1]):
            assert(tools.correlation(imputed[:,i], self.module.site_imputation[:,i]))
        
        for i in range(imputed.shape[0]):
            assert(tools.correlation(imputed[i,:], self.module.site_imputation[i,:]))

        logging.info("PASS")



    # def test_impute_missing_samples_and_snps(self):
    #     logging.info("Testing test_impute_missing_samples_and_snps...")
    #     imputed = loadtxt(self.MISSING_SAMPLES_AND_SNPS_MANUALLY_IMPUTED_FILE)
    #     self.module.impute(self.MIN_SCORE, self.SNPS_FILE, self.MISSING_SAMPLES_AND_SNPS_GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        
    #     assert(imputed.shape[0] == self.module.site_imputation.shape[0])
    #     assert(imputed.shape[1] == self.module.site_imputation.shape[1])

    #     for i in range(imputed.shape[1]):
    #         assert(tools.correlation(imputed[:,i], self.module.site_imputation[:,i]))
        
    #     for i in range(imputed.shape[0]):
    #         assert(tools.correlation(imputed[i,:], self.module.site_imputation[i,:]))
    #     logging.info("PASS")
