"""""
ReFactor
"""""
import logging
import os
import sys
import argparse
from modules import refactor

class RefactorParser( object ):

    ALL_ARGS = ['--pheno', '-k', '-t', '--covar', '--numcomp', '--fs']

    @staticmethod
    def init_args( parser ):
        refactor = parser.add_argument_group('refactor', 'add refactor description here')

        refactor.add_argument('--pheno',    type = str,  help = "a phenotype file for the association test (can include multiple phenotypes)")
        refactor.add_argument('-k',         type = int,   help = "the number of assumed cell types (must be given if --refactor is used)")
        refactor.add_argument('-t',         type = int, default = 500,  help = "the number of sites refactor uses for computing the components. Default is 500") #   monte_carlo_size
   
        #TODO what type is this
        refactor.add_argument('--covar',                                help="for including covariates; these will be used for computing the refactor components as well as for running association test (if --linreg or --logreg are selected)")
        refactor.add_argument('--numcomp',  type = int,                 help = "the number of refactor components to output. Default is K")

        refactor.add_argument('--fs',       type = str, default = 'normal',  help = "feature selection mode (for the first step of refactor)  - normal (defualt) (the standard refactor feature selection; default value), controls (by using the controls only; possible only if a binary phenotype is provided), phenotype (by using a subspace orthogonal to the phenotype)")


    def __init__(self, args, data, output_perfix = ""):
        # validate required args
        if not args.k:
            logging.error("-k option must be specified with the --refactor flag")
            sys.exit(2) 

        # init module to validate other args format
        self.module  = refactor.Refactor(O = data, 
                              K = args.k, 
                              t = args.t, 
                              feature_selection = args.fs.lower().strip(), 
                              num_components = args.numcomp, 
                              phenofile = args.pheno,
                              covar = args.covar,
                              ranked_output_filename = output_perfix + refactor.RANKED_FILENAME, 
                              components_output_filename  = output_perfix + refactor.COMPONENTS_FILENAME)

    def run( self ):
        self.module.run()

