import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from refactor_parser import RefactorParser
from modules import lmm, methylation_data, ewas
from numpy import savetxt, column_stack, array

LMM_OUT_SUFFIX = ".lmm.txt"

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
        lmm_parser.add_argument('--norm', action='store_true', help='Supply this flag in order to normalize covariates matrix (if this flag is not supplied the matrix is not normalized)')
        lmm_parser.add_argument('--calc', action='store_true', help='Supply this float in order to generate logdelta for each site seperatly (if not specified - one logdelta will be generated. Note that this flag has no meaning when --logdelta is supplied too')
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

        if args.logdelta and args.calc:
          logging.warning("LMM will use the same logdelta for all sites. logdelta that was supplied is %s" % args.logdelta)

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

            if args.calc: # run lmm for each site so logdelta will be calculated for each site (TODO move this option as an argument of LMM class)
              cpgnames=[]
              pvalues=[]
              intercepts_betas=[]
              covars_betas=[]
              sites_betas=[]
              sigmas_g=[]
              sigmas_e=[]
              stats=[]

              for i in range(meth_data.sites_size):
                data_site_i = data[:,i].reshape((-1,1)) # n samples by 1 site
                res = module.run(data_site_i, pheno, covars, [meth_data.cpgnames[i]], args.norm, args.logdelta, args.reml)
                cpgname, pvalue, intercept_beta, covariates_beta, site_beta, sigma_e, sigma_g, statistic = res
                
                cpgnames.append(cpgname[0])
                pvalues.append(pvalue[0])
                intercepts_betas.append(intercept_beta[0])
                covars_betas.append(covariates_beta[0])
                sites_betas.append(site_beta[0])
                sigmas_e.append(sigma_e[0])
                sigmas_g.append(sigma_g[0])
                stats.append(statistic[0])

            else: # run lmm on all data - logdelta is calculated once.
              #run lmm
              lmm_results = module.run(data, pheno, covars, meth_data.cpgnames, args.norm, args.logdelta, args.reml)
              cpgnames, pvalues, intercepts_betas, covars_betas, sites_betas, sigmas_e, sigmas_g, stats =  lmm_results
            

            # generate result - by EWAS output format
            ewas_res = ewas.EWASResultsCreator("LMM", array(cpgnames), array(pvalues), statistic = array(stats),\
                                              intercept_coefs = array(intercepts_betas), covars_coefs = array(covars_betas), \
                                              site_coefs = array(sites_betas), sigma_g = array(sigmas_g), sigma_e = array(sigmas_e))

            # save results
            output_file = "results" + LMM_OUT_SUFFIX if output_perfix is None else output_perfix + LMM_OUT_SUFFIX
            ewas_res.save(output_file)

            return ewas_res


        except Exception :
            logging.exception("in lmm")
            raise
