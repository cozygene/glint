import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from refactor_parser import RefactorParser
from modules import lmm, methylation_data
from numpy import savetxt, column_stack

LMM_OUT_SUFFIX = "lmm.txt"

class LMMParser(ModuleParser):
    def __init__(self, parser):
        self.parser = parser
        lmm_parser = parser.add_argument_group('lmm', 'TODO Elior,add lmm description here')

        def kinship_value(val):
            if val not in lmm.AVAILABLE_KINSHIPS:
                common.terminate("kinship must be \"%s\"" % "\" or \"".join(lmm.AVAILABLE_KINSHIPS))
            return val
        lmm_parser.add_argument('--kinship', type = kinship_value, required = True, help = "The way to generate kinship matrix. Options are %s. if \"refactor\" is selected than refactor flags need to be supplied (more info udser \"refactor\" section). " % str(lmm.AVAILABLE_KINSHIPS))
        
        def reml_value(val):
            val = int(val)
            if val != 0 and val != 1:
                common.terminate("reml must be 0 or 1")
            return val
        lmm_parser.add_argument('--reml', type=reml_value, default=1, help='type 1 to use REML (restricted maximum likelihood) or 0 to use ML. Default is 1 (REML)')
        lmm_parser.add_argument('--logdelta', type=float, default=None, help='The value of log(delta) to use. Infers it from the data by default (if not specified)')
        
        # def pc_num_value(val):
        #     val = int(val)
        #     if val < 0:
        #         common.terminate("numpccovar must positive integer (0+)")
        #     return val
        # todo use --numcomp??
        # lmm_parser.add_argument('--numpccovar', metavar='numPCCovars', type=pc_num_value, default=0, help='The number of principal components to use as covariates. If 0 not used as covariates. Default is 0')
        
        super(LMMParser, self).__init__(lmm_parser)
        

    def validate_args(self, args):
        # argument pheno is required for all ewas tests - it can be supplied through --pheno flag of .glint meth data file
        # So, if the datafile supplied is not .glint file - pheno must be supplied as a flag 
        if not args.datafile.name.endswith(methylation_data.GLINT_FORMATTED_EXTENSION):
            self.required_args.append('pheno')
        
        if args.kinship == 'refactor':
            self.refactor = RefactorParser(self.parser)
            self.required_args.extend(self.refactor.required_args)
            self.all_args.extend(self.refactor.all_args)

        super(LMMParser, self).validate_args(args)
      
    def run(self, args, meth_data, output_perfix):
        try:
            if meth_data.phenotype is None and args.pheno is None:
                common.terminate("phenotype file wasn't supplied")

            kinship_data = None

            if args.kinship == 'refactor': # kinship and data to test are the same
                # todo if --lmm provoded woth --refactor there is no need to run refactor twice in order to find ranked sites.
                logging.info("Running lmm with refactor kinship...")
                refactor_meth_data = meth_data.copy() #todo need to copy?
                self.refactor.run(args, refactor_meth_data, output_perfix)

                logging.info("using best %s sites suggested by refactor as data for kinship..." % args.t)
                t_best_sites = self.refactor.module.ranked_sites[:args.t]

                data_for_kinship = meth_data.copy() #todo need to copy?
                data_for_kinship.include(t_best_sites)
                
                # todo handle if sample not in phenotype

                # all data is of dimensions n samplesX m sites
                kinship_data = data_for_kinship.data.transpose()
                kinship = lmm.KinshipCreator(kinship_data, is_normalized = False).create_standard_kinship()

            # all data is of dimensions n samplesX m sites
            covars = meth_data.covar # no need to transpose since covars are nXm (n samples)
            data = meth_data.data.transpose() # data to test
            pheno = meth_data.phenotype #should transpose? todo

            # initialize lmm with kinship
            module = lmm.LMM(kinship)
            sorted_cpgnames, sorted_cpg_indices, pvalues = module.run(data, pheno, covars, meth_data.cpgnames, False, args.logdelta, args.reml)

            output_file = LMM_OUT_SUFFIX if output_perfix is None else output_perfix + LMM_OUT_SUFFIX
            logging.info("saving LMM output to file %s" % output_file)
            savetxt(output_file, column_stack((sorted_cpgnames, pvalues)), fmt = '%s16')

        except Exception :
            logging.exception("in lmm")
            raise
