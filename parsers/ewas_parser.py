import sys
import logging
from utils import common
"""""
EWAS
"""""
class EWASParser(object):
    ALL_ARGS = ['--pheno', '--covar', '--linreg', '--logreg']
   
    @staticmethod
    def init_args( parser ):
        all_args = []
        ewas = parser.add_argument_group('ewas', 'add ewas description here')
        ewas.add_argument('--pheno', type = str, help = "A phenotype file")
        ewas.add_argument('--covar', type = str, help="A covariates file")

        group_reg = ewas.add_mutually_exclusive_group(required = False) #TODO this is not working!
        group_reg.add_argument('--linreg', help = "Run a linear regression analysis; --pheno must be provided (executed by default if --ewas is selected)")
        group_reg.add_argument('--logreg', help = "Run a logistic regression analysis; --pheno must be provided and be a binary phenotype")
        
    def __init__( self, args ):

        ## put here parser checks
        ## call here ewas so it check its own


        if args.logreg and not args.pheno: # this won't happen since pheno is required
            logging.error("must specify phenofile with logreg")
            common.terminate(self.__class__.__name__)

        if args.linreg and not args.pheno:
            logging.error("must specify phenofile with linreg")
            common.terminate(self.__class__.__name__)

        self.module = None

    def run(self):
        # self.module.run()
        pass