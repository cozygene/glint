import sys
import logging
"""""
EWAS
"""""
class EWASParser(object):
    ALL_ARGS = ['--pheno', '--covar', '--linreg', '--logreg']
   
    @staticmethod
    def init_args( parser ):
        all_args = []
        ewas = parser.add_argument_group('ewas', 'add ewas description here')
        ewas.add_argument('--pheno',     type = str,       help = "skkk")
        ewas.add_argument('--covar',     type = str,        help="for including covariates; these will be used for computing the refactor components as well as for running association test (if --linreg or --logreg are selected)")


        group_reg = ewas.add_mutually_exclusive_group(required = False) #TODO this is not working!
        group_reg.add_argument('--linreg',        help = "performs a linear regression. If selecting that option must provide a phenotype (--phenofile)")
        group_reg.add_argument('--logreg',        help = "performs a logistic regression. If selecting that option must get a phenotype (--phenofile), and it must be a binary phenotype.")
        
    def __init__( self, args ):

        ## put here parser checks
        ## call here ewas so it check its own


        if args.logreg and not args.pheno: # this won't happen since pheno is required
            logging.error("must specify phenofile with logreg")
            sys.exit(2)

        if args.linreg and not args.pheno:
            logging.error("must specify phenofile with linreg")
            sys.exit(2)

        self.module = None

    def run(self):
        # self.module.run()
        pass