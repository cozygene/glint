import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import predictor

                

SITES_SCORES_FILE = 'parsers/assets/sites_scores_list'
SITES_SNPS_FILE = 'parsers/assets/site_snps_list'
SITES_IDS_FILE ='parsers/assets/sites_ids_list'
SNPS_IDS_FILE = 'parsers/assets/snps_ids_list'
SITES_SNPS_COEFF_FILE = 'parsers/assets/sites_snps_coeff_list'

class PredictorParser(ModuleParser):
    def __init__(self, parser):

        predictor = parser.add_argument_group('predictor', 'predict methylation levels by SNPs. TODO Elior,add ewas description here')
        predictor.add_argument('--psnp',  type = argparse.FileType('r'), required = True, help = "plink .snp file")
        predictor.add_argument('--pgeno',  type = argparse.FileType('r'), required = True, help = "pling .geno file")
        predictor.add_argument('--pind',  type = argparse.FileType('r'), required = True, help = "plink .ind file")
        predictor.add_argument('--score',  type = float, default = 0.5, help = "a score (between 0 and 1) - predict only cpgs with prediction correlation score at least S ")
        predictor.add_argument('--maxmiss',  type = float, default = 0.03, help = "The maximal amount of missing values allowed in the plink data files (percentage, number between 0 and 1). if sample X has more than <minmiss> percentage missing values (out of all it's snps) dont predict it, the same for snp")
        super(PredictorParser, self).__init__(predictor)
        

    def run(self, args):
        try:
            self.module  = predictor.Predictor(SITES_SCORES_FILE, SITES_SNPS_FILE, SITES_IDS_FILE, SNPS_IDS_FILE, SITES_SNPS_COEFF_FILE)
            self.module.predict(args.score, args.psnp, args.pgeno, args.pind, args.maxmiss)
            meth_data = self.module.meth_data()
            meth_data.save('predicted')
        except Exception :
            logging.exception("in predictor")
            raise