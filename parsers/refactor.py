"""""
ReFactor
"""""
import logging
import os
import sys
from modules import refactor

def init_parser( parser ):
    refactor = parser.add_argument_group('refactor', 'add refactor description here')
    refactor.add_argument('-k',         type = int,                 help = "the number of assumed cell types (must be given if --refactor is used)") 
    refactor.add_argument('-t',         type = int, default = 500,  help = "the number of sites refactor uses for computing the components. Default is 500") #   monte_carlo_size
    #TODO what type is this
    refactor.add_argument('--covar',                                help="for including covariates; these will be used for computing the refactor components as well as for running association test (if --linreg or --logreg are selected)")
    refactor.add_argument('--numcomp',  type = int,                 help = "the number of refactor components to output. Default is K")
    refactor.add_argument('--pheno',    type = str,                 help = "a phenotype file for the association test (can include multiple phenotypes)")
    # TODO can be together linreg + logreg
    refactor.add_argument('--linreg',   action='store_true',        help = "performs a linear regression. If selecting that option must provide a phenotype (--phenofile)")
    refactor.add_argument('--logreg',    action='store_true',        help = "performs a logistic regression. If selecting that option must get a phenotype (--phenofile), and it must be a binary phenotype.")
    refactor.add_argument('--fs',       type = str, default = 'normal',  help = "feature selection mode (for the first step of refactor)  - normal (defualt) (the standard refactor feature selection; default value), controls (by using the controls only; possible only if a binary phenotype is provided), phenotype (by using a subspace orthogonal to the phenotype)")



def parse(data, args):
    if not args.k:
        logging.error("-k option must be specified with the --refactor flag")
        sys.exit(2) 
    
    if not args.numcomp:
        args.numcomp = args.k
    
    if args.pheno and not os.path.exists(args.pheno) :
        logging.error("The file '%s' doesn't exist. Exiting" % args.pheno)
        sys.exit(2)

    if args.logreg and not args.pheno:
        logging.error("must specify phenofile with logreg")
        sys.exit(2)

    if args.logreg and not args.pheno:
        logging.error("must specify phenofile with linreg")
        sys.exit(2)

    # TODO should have default value
    if not args.fs:
        raise Exception("fs")
    elif args.fs == 'normal':
        pass
    elif args.fs == 'phenotype':
        pass
    elif args.fs == 'controls':
        pass
    else:
        logging.error("choose - normal (defualt) (the standard refactor feature selection; default value), controls (by using the controls only; possible only if a binary phenotype is provided), phenotype (by using a subspace orthogonal to the phenotype)")
        sys.exit(2)
    logging.info("GREAT")
    refactor.refactor( O = data, K = args.k, t = args.t, num_components = args.numcomp )
