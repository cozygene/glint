from modules import predictor,methylation_data
from numpy import loadtxt, array_equal
import logging
from tests.test_tools import tools

class PredictorTester():

    SNPS_FILE = "tests/predictor/files/reut22_5snps.snps"
    GENO_FILE = "tests/predictor/files/reut22_10samples_5snps.geno"
    MISSING_SAMPLES_GENO_FILE = "tests/predictor/files/reut22_10samples_5snps-missing_samples.geno"
    MISSING_SNPS_GENO_FILE = "tests/predictor/files/reut22_10samples_5snps-missing_snps.geno"
    MISSING_SAMPLES_AND_SNPS_GENO_FILE = "tests/predictor/files/reut22_10samples_5snps-missing_samples_and_snps.geno"
    IND_FILE = "tests/predictor/files/reut22_10s.ind"
    MANUALLY_PREDICTED_FILE = "tests/predictor/files/reut22_test_out" 
    MISSING_SAMPLES_MANUALLY_PREDICTED_FILE = "tests/predictor/files/reut22_test_out-missing_samples" 
    MISSING_SNPS_MANUALLY_PREDICTED_FILE = "tests/predictor/files/reut22_test_out-missing_snps" 
    MISSING_SAMPLES_AND_SNPS_MANUALLY_PREDICTED_FILE = "tests/predictor/files/reut22_test_out-missing_samples_and_snps"
    SITES_SCORES_FILE = 'parsers/assets/sites_scores_list'
    SITES_SNPS_FILE = 'parsers/assets/site_snps_list'
    SITES_IDS_FILE ='parsers/assets/sites_ids_list'
    SNPS_IDS_FILE = 'parsers/assets/snps_ids_list'
    SITES_SNPS_COEFF_FILE = 'parsers/assets/sites_snps_coeff_list'

    MIN_SCORE = 0.0005
    MIN_MISSING_VALUES = 0.2
    def __init__(self):
        logging.info("Testing Started on PredictorTester")
        self.module  = predictor.Predictor(self.SITES_SCORES_FILE, self.SITES_SNPS_FILE, self.SITES_IDS_FILE, self.SNPS_IDS_FILE, self.SITES_SNPS_COEFF_FILE)
        self.test_predict()
        self.test_predict_missing_samples()
        self.test_predict_missing_snps()
        self.test_predict_missing_samples_and_snps()
        logging.info("Testing Finished on PredictorTester")

    def test_predict(self):
        logging.info("Testing test_predict...")
        predicted = loadtxt(self.MANUALLY_PREDICTED_FILE)
        self.module.predict(self.MIN_SCORE, self.SNPS_FILE, self.GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        
        assert(predicted.shape[0] == self.module.site_prediction.shape[0])
        assert(predicted.shape[1] == self.module.site_prediction.shape[1])
        for i in range(predicted.shape[1]):
            assert(tools.correlation(predicted[:,i], self.module.site_prediction[:,i]))
        
        for i in range(predicted.shape[0]):
            assert(tools.correlation(predicted[i,:], self.module.site_prediction[i,:]))
        logging.info("PASS")

    def test_predict_missing_samples(self):
        logging.info("Testing test_predict_missing_samples...")
        predicted = loadtxt(self.MISSING_SAMPLES_MANUALLY_PREDICTED_FILE)
        self.module.predict(self.MIN_SCORE, self.SNPS_FILE, self.MISSING_SAMPLES_GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        assert(predicted.shape[0] == self.module.site_prediction.shape[0])
        assert(predicted.shape[1] == self.module.site_prediction.shape[1])
        for i in range(predicted.shape[1]):
            assert(tools.correlation(predicted[:,i], self.module.site_prediction[:,i]))
        
        for i in range(predicted.shape[0]):
            assert(tools.correlation(predicted[i,:], self.module.site_prediction[i,:]))
        logging.info("PASS")


    def test_predict_missing_snps(self):
        logging.info("Testing test_predict_missing_snps...")
        predicted = loadtxt(self.MISSING_SNPS_MANUALLY_PREDICTED_FILE)
        self.module.predict(self.MIN_SCORE, self.SNPS_FILE, self.MISSING_SNPS_GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        
        assert(predicted.shape[0] == self.module.site_prediction.shape[0])
        assert(predicted.shape[1] == self.module.site_prediction.shape[1])

        for i in range(predicted.shape[1]):
            assert(tools.correlation(predicted[:,i], self.module.site_prediction[:,i]))
        
        for i in range(predicted.shape[0]):
            assert(tools.correlation(predicted[i,:], self.module.site_prediction[i,:]))

        logging.info("PASS")



    def test_predict_missing_samples_and_snps(self):
        logging.info("Testing test_predict_missing_samples_and_snps...")
        predicted = loadtxt(self.MISSING_SAMPLES_AND_SNPS_MANUALLY_PREDICTED_FILE)
        self.module.predict(self.MIN_SCORE, self.SNPS_FILE, self.MISSING_SAMPLES_AND_SNPS_GENO_FILE, self.IND_FILE, self.MIN_MISSING_VALUES)
        
        assert(predicted.shape[0] == self.module.site_prediction.shape[0])
        assert(predicted.shape[1] == self.module.site_prediction.shape[1])

        for i in range(predicted.shape[1]):
            assert(tools.correlation(predicted[:,i], self.module.site_prediction[:,i]))
        
        for i in range(predicted.shape[0]):
            assert(tools.correlation(predicted[i,:], self.module.site_prediction[i,:]))
        logging.info("PASS")
