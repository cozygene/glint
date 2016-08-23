import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from lmm_parser import LMMParser
from modules import ewas, methylation_data



LINREG_OUT_SUFFIX = ".glint.linreg.txt"
LOGREG_OUT_SUFFIX = ".glint.logreg.txt"

"""""
EWAS
"""""
class EWASParser(ModuleParser):
    def __init__(self, parser):
        ewas = parser.add_argument_group('ewas', 'TODO Elior,add ewas description here')

        # Note that argument '--pheno' is required for all EWAS tests. but dont add it to dependencies list (dependencies = ['--pheno'])
        # since it can be supplied through the meth_data object (if .glint file was provided and not meth data matrix)
        ewas.add_argument('--linreg', action = "store_true", help = "Run a linear regression analysis (executed by default if --ewas is selected)")
        ewas.add_argument('--logreg', action = "store_true", help = "Run a logistic regression analysis")
        # Note: lmm module is handled not the very best way since there was no time. it appears here under "EWAS" but the glin.py handles it as 
        # an independed module
        ewas.add_argument('--lmm', dependencies = ["--ewas"], action = "store_true", help = "Run a linear mixed model test. More explanation and options described under \"lmm\"")
        
        self.lmm_parser = LMMParser(parser)
        super(EWASParser, self).__init__(ewas)

    def validate_args(self, args):
        # argument pheno is required for all ewas tests - it can be supplied through --pheno flag of .glint meth data file
        # So, if the datafile supplied is not .glint file - pheno must be supplied as a flag 
        if not args.datafile.name.endswith(methylation_data.GLINT_FORMATTED_EXTENSION):
            self.required_args.append('pheno')
        
        # make sure user choose only one EWAS test (mutually exlusive group is not supported...)
        test_counter = 0
        if args.lmm:
            test_counter += 1
        if args.logreg:
            test_counter += 1
        if args.linreg:
            test_counter += 1
        if args.wilc:
            test_counter +=1
        if test_counter > 1:
            common.terminate("Choose only one EWAS test.")

        # add lmm parser if lmm was chosen
        if args.lmm:
            self.lmm_parser.validate_args(args)
            self.all_args.extend(self.lmm_parser.all_args)

        # default test is linear regression
        if test_counter == 0:
            args.linreg = True
            logging.info("No EWAS test was chosen, running linerar regression by default.")


        super(EWASParser, self).validate_args(args)

    def runLMM(self, args, meth_data):
        return self.lmm_parser.run(args = args,
                                   meth_data = meth_data,
                                   output_perfix = args.out)
    
    def runRegression(self, meth_data, regression_class, test_name, output_file):
        module = regression_class(methylation_data = meth_data)

        #run lin reg and save result by EWAS output format
        results = module.run()
        cpgnames, pvalues, fstats, intercept_beta, covars_betas, site_beta = results

        ewas_res = ewas.EWASResultsCreator(test_name, cpgnames, pvalues, statistic = fstats,              \
                                          intercept_coefs = intercept_beta, covars_coefs = covars_betas, \
                                          site_coefs = site_beta)

        # save results
        ewas_res.save(output_file)
        return ewas_res

    def runLinReg(self, args, meth_data):
        output_perfix = args.out
        output_file = "results" + LINREG_OUT_SUFFIX if output_perfix is None else output_perfix + LINREG_OUT_SUFFIX
        return self.runRegression(meth_data, ewas.LinearRegression, "LinReg", output_file)

    def runLogReg(self, args, meth_data):
        output_perfix = args.out
        output_file = "results" + LOGREG_OUT_SUFFIX if output_perfix is None else output_perfix + LOGREG_OUT_SUFFIX
        return self.runRegression(meth_data, ewas.LogisticRegression, "LogReg", output_file)

    def run(self, args, meth_data):
        try:
            if meth_data.phenotype is None and args.pheno is None:
                common.terminate("phenotype file wasn't supplied")

            # ewas test must be called after refactor
            if args.lmm:
                return self.runLMM(args, meth_data)
                
            if args.linreg:
                return self.runLinReg(args, meth_data)

            if args.logreg:
                return self.runLogReg(args, meth_data)

        except Exception :
            logging.exception("in ewas")
            raise

