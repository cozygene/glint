import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import ewas
"""""
EWAS
"""""
class EWASParser(ModuleParser):
    def __init__(self, parser):

        ewas = parser.add_argument_group('ewas', 'TODO Elior,add ewas description here')
        ewas.add_argument('--linreg', action='store_const', const='linear_regression',   help = "Run a linear regression analysis; --pheno must be provided (executed by default if --ewas is selected)")
        ewas.add_argument('--logreg', action='store_const', const='logistic_regression', help = "Run a logistic regression analysis; --pheno must be provided and be a binary phenotype")
        
        super(EWASParser, self).__init__(ewas)
        

    def validate_args(self, args):
        super(EWASParser, self).validate_args(args)
        # default test is linear regression
        if not args.logreg and not args.linreg:
            self.tests = ['linear_regression']
        else:
            self.tests = []
            if args.logreg is not None:
                self.tests.append(args.logreg)
            if args.linreg is not None:
                self.tests.append(args.linreg)


    def run(self, args, meth_data, output_perfix = ""):
        try:
            if meth_data.phenotype is None and args.pheno is None:
                common.terminate("phenotype file wasn't supplied")

            self.module  = ewas.EWAS(methylation_data = meth_data, tests_list = self.tests)
            self.module.run()
        except Exception :
            logging.exception("in ewas")
            raise