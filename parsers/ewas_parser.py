import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from lmm_parser import LMMParser
from modules import ewas, methylation_data


"""""
EWAS
"""""
class EWASParser(ModuleParser):
    def __init__(self, parser):
        ewas = parser.add_argument_group('ewas', 'TODO Elior,add ewas description here')

        # Note that argument '--pheno' is required for all EWAS tests. but dont add it to dependencies list (dependencies = ['--pheno'])
        # since it can be supplied through the meth_data object (if .glint file was provided and not meth data matrix)
        ewas.add_argument('--linreg', action='store_const', const='linear_regression',   help = "Run a linear regression analysis (executed by default if --ewas is selected)")
        ewas.add_argument('--logreg', action='store_const', const='logistic_regression', help = "Run a logistic regression analysis")
        # Note: lmm module is handled not the very best way since there was no time. it appears here under "EWAS" but the glin.py handles it as 
        # an independed module
        ewas.add_argument('--lmm', dependencies = ["--ewas"], action='store_const', const='logistic_regression', help = "Run a linear mixed model test. More explanation and options described under \"lmm\"")
        
        self.lmm_parser = LMMParser(parser)
        super(EWASParser, self).__init__(ewas)

    def validate_args(self, args):
        # argument pheno is required for all ewas tests - it can be supplied through --pheno flag of .glint meth data file
        # So, if the datafile supplied is not .glint file - pheno must be supplied as a flag 
        if not args.datafile.name.endswith(methylation_data.GLINT_FORMATTED_EXTENSION):
            self.required_args.append('pheno')

        if args.lmm:
            self.lmm_parser.validate_args(args)
            self.all_args.extend(self.lmm_parser.all_args)
        # default test is linear regression
        if not args.logreg and not args.lmm:
            logging.info("No EWAS test was chosen, running linerar regression by default.")
            self.tests = ['linear_regression']
        else:
            self.tests = []
            if args.logreg is not None:
                self.tests.append(args.logreg)
            if args.linreg is not None:
                self.tests.append(args.linreg)

        super(EWASParser, self).validate_args(args)

    def run(self, args, meth_data):
        try:
            if meth_data.phenotype is None and args.pheno is None:
                common.terminate("phenotype file wasn't supplied")

            # ewas test must be called after refactor
            if args.lmm:
                #return lmm results object
                return self.lmm_parser.run(args = args,
                                           meth_data = meth_data,
                                           output_perfix = args.out)
            else:
                self.module  = ewas.EWAS(methylation_data = meth_data, tests_list = self.tests)
                self.module.run()
        except Exception :
            logging.exception("in ewas")
            raise

