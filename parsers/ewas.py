"""""
EWAS
"""""
def init_parser( parser ):

    ewas = parser.add_argument_group('ewas', 'add ewas description here')
    ewas.add_argument('--covar',                                help="for including covariates; these will be used for computing the refactor components as well as for running association test (if --linreg or --logreg are selected)")
    ewas.add_argument('--pheno',    type = str,                 help = "a phenotype file for the association test (can include multiple phenotypes)")
    
def parse(args):
    pass