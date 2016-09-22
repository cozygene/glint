import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from lmm_parser import LMMParser
from modules import ewas, methylation_data



LINREG_OUT_SUFFIX = ".glint.linreg.txt"
LOGREG_OUT_SUFFIX = ".glint.logreg.txt"
WILCOXON_OUT_SUFFIX = ".glint.wilcoxon.txt"

"""""
EWAS
"""""
class EWASParser(ModuleParser):
    def __init__(self, parser):
        """
        ewas allows the following tests (it allows to execute only one test at a time):
         - stdth - as described in the help
         - covar - list of the names of covariates to use.
                if flag is not set will not use any covariates. if set and no name if specified will use all of the covariates. o
                therwise (names specified) will only use the covariates which name is in the list
         - pheno - list of names of phenotypes to use. behaves liske --covar
         -stdth threshold for excluding low variance sites (all sites with std lower than this threshold will be excluded
         - linear regression - the default test
         - logistic regression (validates that phenotype is  binary)
         - wilcoxon rank-sum (validates that phenotype is binary, terminates if --covar flag was supplied (but not if there is covaraites in a glint file))
         - lmm 
        * terminates if phenotype file wasn't supplied (with --pheno of with glint file)

        * output file is different for each test:
            
          -for linear and logistic regression:
            LinReg/LogReg:ID (cpgnames), chromosome, MAPINFO (position), p-value, q-value, intercept (the intercept coefficient), V1 (first covar coefficient),...
            , Vn (last covar coefficient), beta (site under test coefficient), statistic, UCSC_RefGene_Name (gene), Relation_to_UCSC_CpG_Island (category)
           
          - for wilcoxon test
            Wilcoxon:ID (cpgnames), chromosome, MAPINFO (position), p-value, q-value, statistic, UCSC_RefGene_Name (gene), Relation_to_UCSC_CpG_Island (category)
          
          - for lmm see LMMParser documentation
        
        * plot 
            in order to plot the output call --plot with the plot you want.
            you can also execute plots after the test by supplying the test's result file
        """
        ewas = parser.add_argument_group('ewas', 'TODO Elior,add ewas description here')

        # phenotype is required for all EWAS tests
        ewas.add_argument('--pheno', required = True, type = str, nargs='*', help = "list of phenotypes names to use. If no name is specified will use all the phenotypes. If flag is not set, will not use any phenotype")
        # covar is for lmm, linreg and logreg tests
        ewas.add_argument('--covar', type = str, nargs='*', help = "list of covariates names to use. If no name is specified will use all the covariates. If flag is not set, will not use any covariate")
        def std_value(num):
            try:
                num = float(num)
            except:
                common.terminate("minstd must be a float between 0 and 1")
            if not (num <= 1 and num >=0):
                common.terminate("minstd must be a float between 0 and 1")
            return num
        ewas.add_argument('--stdth',  type = std_value, help = "threshold for excluding low variance sites (all sites with std lower than this threshold will be excluded)")  
      
        ewas.add_argument('--linreg', action = "store_true", help = "Run a linear regression analysis (executed by default if --ewas is selected)")
        ewas.add_argument('--logreg', action = "store_true", help = "Run a logistic regression analysis")
        ewas.add_argument('--wilc',   action = "store_true", help = "Run Wilcoxon rank-sum test")
        # Note: lmm module is handled not the very best way since there was no time. it appears here under "EWAS" but the glin.py handles it as 
        # an independed module
        ewas.add_argument('--lmm', dependencies = ["--ewas"], action = "store_true", help = "Run a linear mixed model test. More explanation and options described under \"lmm\"")
        
        self.lmm_parser = LMMParser(parser)
        super(EWASParser, self).__init__(ewas)

    def validate_args(self, args):
        # make sure user choose only one EWAS test (mutually exlusive group is not supported...)
        # argument pheno is required for all ewas tests - it can be supplied through --pheno flag of .glint meth data file
        # So, if the datafile supplied is not .glint file - pheno must be supplied as a flag         
        if not args.datafile.name.endswith(methylation_data.GLINT_FORMATTED_EXTENSION):
            self.required_args.append('phenofiles')

        super(EWASParser, self).validate_args(args)
        
        if len(args.pheno) > 1: # 0 is all phenotypes in the data (which could be one)
            common.terminate("must supply only one phenotype for EWAS")

        test_counter = 0
        if args.lmm:
            test_counter += 1
        if args.logreg:
            test_counter += 1
        if args.linreg:
            test_counter += 1
        if args.wilc:
            if args.covar:
                common.terminate("Wilcoxon test cannot take any covaraites. remove --covar flag")
            test_counter +=1
        if test_counter > 1:
            common.terminate("Choose only one EWAS test.")

        # add lmm parser if lmm was chosen
        if args.lmm:
            self.lmm_parser.validate_args(args)
            self.all_args.extend(self.lmm_parser.all_args)
            self.required_args.extend(self.lmm_parser.required_args)

        # default test is linear regression
        if test_counter == 0:
            args.linreg = True
            logging.info("No EWAS test was chosen, running linerar regression by default.")


        

    def runLMM(self, args, meth_data, pheno, covars):
        return self.lmm_parser.run(args = args,
                                   meth_data = meth_data,
                                   output_perfix = args.out,
                                   covars = covars,
                                   pheno = pheno)

    def runRegression(self, data, regression_class, test_name, output_file, cpgnames, pheno, covars = None):
        module = regression_class(data, cpgnames, pheno, covars)

        #run lin reg and save result by EWAS output format
        results = module.run()
        cpgnames, pvalues, fstats, intercept_beta, covars_betas, site_beta = results

        ewas_res = ewas.EWASResultsCreator(test_name, cpgnames, pvalues, statistic = fstats,              \
                                          intercept_coefs = intercept_beta, covars_coefs = covars_betas, \
                                          site_coefs = site_beta)

        # save results
        ewas_res.save(output_file)
        return ewas_res

    def runLinReg(self, args, data, cpgnames, pheno, covars):
        output_perfix = args.out
        output_file = "results" + LINREG_OUT_SUFFIX if output_perfix is None else output_perfix + LINREG_OUT_SUFFIX
        return self.runRegression(data, ewas.LinearRegression, "LinReg", output_file, cpgnames, pheno, covars)

    def runLogReg(self, args, data, cpgnames, pheno, covars):
        output_perfix = args.out
        output_file = "results" + LOGREG_OUT_SUFFIX if output_perfix is None else output_perfix + LOGREG_OUT_SUFFIX
        return self.runRegression(data, ewas.LogisticRegression, "LogReg", output_file, cpgnames, pheno, covars)

    def runWilcoxon(self, args, data, cpgnames, pheno):
        output_perfix = args.out
        output_file = "results" + WILCOXON_OUT_SUFFIX if output_perfix is None else output_perfix + WILCOXON_OUT_SUFFIX
        test_module = ewas.Wilcoxon(data, cpgnames, pheno)
        cpgnames, pvalues, fstats = test_module.run()
        ewas_res = ewas.EWASResultsCreator("Wilcoxon", cpgnames, pvalues, statistic = fstats)
        ewas_res.save(output_file)
        return ewas_res


    def run(self, args, meth_data):
        try:
            if args.stdth is not None:
                meth_data.remove_lowest_std_sites(args.stdth)

            pheno = meth_data.get_phenotype_subset(args.pheno)
            if (pheno.shape[1] != 1): # check if selected more than one phenotype
                common.terminate("must supply only one phenotype for EWAS")

            covars = meth_data.get_covariates_subset(args.covar)

            # ewas test must be called after refactor
            if args.lmm:
                return self.runLMM(args, meth_data, pheno, covars)
                
            if args.linreg:
                return self.runLinReg(args, meth_data.data, meth_data.cpgnames, pheno, covars)

            if args.logreg:
                return self.runLogReg(args, meth_data.data, meth_data.cpgnames, pheno, covars)

            if args.wilc:
                return self.runWilcoxon(args, meth_data.data, meth_data.cpgnames, pheno)

        except Exception :
            logging.exception("in ewas")
            raise

