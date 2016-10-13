import os
import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import imputing

SITES_SCORES_FILE = os.path.join(os.path.dirname(__file__), 'assets/sites_scores_list')
SITES_SNPS_FILE = os.path.join(os.path.dirname(__file__), 'assets/site_snps_list')
SITES_IDS_FILE = os.path.join(os.path.dirname(__file__), 'assets/sites_ids_list')
SNPS_IDS_FILE = os.path.join(os.path.dirname(__file__), 'assets/snps_ids_list')
SITES_SNPS_COEFF_FILE = os.path.join(os.path.dirname(__file__), 'assets/sites_snps_coeff_list')

class ImputingParser(ModuleParser):
    def __init__(self, parser):
        """
        Notes:
        impute methylation level of methylation sites by snps.
    
        - input: three EIGENSTRAT files not (seperated by chromosome!):
            - .snp file 
            - .geno file
            - .ind file
        - output: imputed methylationData object

        - ignores G/C and T/A SNPs (there is no flag to disable this right now)
        - ignores samples with more than --maxmiss snps that are missing (missing values).
        - ignores SNPs with more than --maxmiss samples that are missing (missing values)
        - ignores sites with low score, uses only sites with correlation score at least --score
          (sites scores list is in the file "sites_scores_list")
        - terminates and warns if all sites are ignored.
        - terminates and warns if all samples are ignored.
        - imputes any site S if this site is relevant (wasn't ignored) and if we have information of at least 1 of the SNPs that impute site S.

        - the model used for imputing is encoded in the "site_snps_list" and "sites_snps_coeff_list" files as follows:
          - "site_snps_list" file (the SNPs that imputes the sites):
              the list of numbers appear in the i'th line are list of the SNPs ids* that "explain" the site which id** is i.
              *SNPs ids list is found in "snps_ids_list" file: the name of the SNP which id is j appears at the j'th line in that file
              **sites ids list is found in "sites_ids_list" file: the name of the site which id is i appears at the i'th line in that file
          - "sites_snps_coeff_list" file (the coefficients of the SNPs):
              the numbers at the i'th line are the coefficients of the SNPs which impute ("explain") the i'th site.
              the j'th number in the i'th line is the coefficient of the j'th snp at the i'th line in the file "site_snps_list"
        """
        impute = parser.add_argument_group('imputation', 'Impute methylation levels by SNPs.')
        impute.add_argument('--snp',  type = argparse.FileType('r'), required = True, help = "EIGENSTRAT .snp file")
        impute.add_argument('--geno',  type = argparse.FileType('r'), required = True, help = "EIGENSTRAT .geno file")
        impute.add_argument('--ind',  type = argparse.FileType('r'), required = True, help = "EIGENSTRAT .ind file")
        impute.add_argument('--score',  type = float, default = 0.5, help = "a score (between 0 and 1) - imputes only sites with correlation score at least S")
        impute.add_argument('--maxmiss',  type = float, default = 0.03, help = "The maximal amount of missing values allowed for SNPs in the EIGENSTRAT data files (number between 0 and 1).")
        super(ImputingParser, self).__init__(impute)
        

    def run(self, args, output_perfix = None):
        try:
            self.module  = imputing.Imputation(SITES_SCORES_FILE, SITES_SNPS_FILE, SITES_IDS_FILE, SNPS_IDS_FILE, SITES_SNPS_COEFF_FILE)
            self.module.impute(args.score, args.snp, args.geno, args.ind, args.maxmiss)
            meth_data = self.module.meth_data()
            
            prefix = "imputation"
            if output_perfix is not None:
                prefix = output_perfix + ".imputation"

            meth_data.save_serialized_data(prefix)
        except Exception :
            logging.exception("in imputing")
            raise
