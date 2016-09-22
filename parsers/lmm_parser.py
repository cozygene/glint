import os
import time
import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from refactor_parser import RefactorParser
from modules import lmm, methylation_data, ewas
from numpy import savetxt, column_stack, array, loadtxt

LMM_OUT_SUFFIX = ".glint.lmm.txt"

class LMMParser(ModuleParser):
    def __init__(self, parser):
        """
        LMM test
        --kinship:  today there are two onptions for a kinship:
                        - a path to a file with the kinship matrix
                        - "refactor" to run refactor and use the suggested best t sites as a kinship
                           * Note that if user selected refactor, than all refactor flags will be avaliable for him.
                             if there is a flag required for refactor - it will be required for the kinship too (today, only k)
                           * Note that if the user run both refactor and lmm with refactor kinship 
                                glint.py --refactor --datafile <> --ewas --lmm --kinship refactor
                             than refactor will be executed twice.
        --reml:     whether to use REML (restricted maximum likelihood)  or ML (maximum likelihood)
                    1 is for REML (default)
                    0 is for ML
        --norm:     if this flag is selectes than the covariates will be normalized before running LMM (if this flag is not supplied the matrix is not normalized)
        --oneld     select this flag to generate one logdelta value for all the sites (by default, without hte flag, it is generated seperatly for each site)
        * terminates if phenotype file wasn't supplied (with --pheno of with glint file)

        * output file is a matrix with the following columns ( at this order as for today ):
            LMM:ID (cpgnames), chromosome, MAPINFO (position), p-value, q-value, intercept (the intercept coefficient), V1 (first covar coefficient),...
            , Vn (last covar coefficient), beta (site under test coefficient), statistic, sigma-e, sigma-g, UCSC_RefGene_Name (gene), Relation_to_UCSC_CpG_Island (category)
        
        * plot 
            in order to plot the output call --plot with the plot you want.
            you can also execute plots after the test by supplying the test's result file
        """
        self.parser = parser
        lmm_parser = parser.add_argument_group('lmm', 'TODO Elior,add lmm description here')

        def kinship_value(val):
            if val not in lmm.AVAILABLE_KINSHIPS and not os.path.exists(val):
                common.terminate("kinship must be a matrix file or  \"%s\"" % "\" or \"".join(lmm.AVAILABLE_KINSHIPS))
            if os.path.exists(val):
                val = open(val, 'rb')
            return val
        lmm_parser.add_argument('--kinship', type = kinship_value, required = True, help = "The way to generate kinship matrix. Options are %s. if \"refactor\" is selected than refactor flags need to be supplied (more info under \"refactor\" section). " % str(lmm.AVAILABLE_KINSHIPS))
        
        def reml_value(val):
            val = int(val)
            if val != 0 and val != 1:
                common.terminate("reml must be 0 or 1")
            return val
        lmm_parser.add_argument('--reml', type=reml_value, default=1, help='type 1 to use REML (restricted maximum likelihood) or 0 to use ML. Default is 1 (REML)')
        lmm_parser.add_argument('--norm', action='store_true', help='Supply this flag in order to normalize covariates matrix (if this flag is not supplied the matrix is not normalized)')
        lmm_parser.add_argument('--oneld', action='store_true', help='select this in order to generate logdelta once for all sites (by default logdelta is generated seperatly for each site)')
        # def pc_num_value(val):
        #     val = int(val)
        #     if val < 0:
        #         common.terminate("numpccovar must positive integer (0+)")
        #     return val
        # should  use --numcomp??
        # lmm_parser.add_argument('--numpccovar', metavar='numPCCovars', type=pc_num_value, default=0, help='The number of principal components to use as covariates. If 0 not used as covariates. Default is 0')
        
        super(LMMParser, self).__init__(lmm_parser)
        

    def validate_args(self, args):
        # argument pheno is required for all ewas tests - it can be supplied through --pheno flag of .glint meth data file
        # So, if the datafile supplied is not .glint file - pheno must be supplied as a flag 
        if not args.datafile.name.endswith(methylation_data.GLINT_FORMATTED_EXTENSION):
            self.required_args.append('phenofiles')
        
        if args.kinship == 'refactor':
            self.refactor = RefactorParser(self.parser)
            self.all_args.extend(self.refactor.all_args)
            self.required_args.extend(self.refactor.required_args)

        super(LMMParser, self).validate_args(args)
      
    def run(self, args, meth_data, pheno, output_perfix, covars = None):
        try:
            kinship_data = None

            if type(args.kinship) == file: #kinship is provided via file
                logging.info("loading kinship from %s" % args.kinship.name)
                kinship = loadtxt(args.kinship)

            elif args.kinship == 'refactor': # kinship and data to test are the same
                # todo if --lmm provided with --refactor there is no need to run refactor twice in order to find ranked sites.
                logging.info("Running lmm with refactor kinship...")
                refactor_meth_data = meth_data.copy()
                self.refactor.run(args, refactor_meth_data, output_perfix)

                logging.info("using best %s sites suggested by refactor as data for kinship..." % args.t)
                t_best_sites = self.refactor.module.ranked_sites[:args.t]
                
                data_for_kinship = meth_data.copy()
                data_for_kinship.include(t_best_sites)
                
                # all data is of dimensions n samplesX m sites
                kinship_data = data_for_kinship.data.transpose()
                kinship = lmm.KinshipCreator(kinship_data, is_normalized = False).create_standard_kinship()

            
                    
            # all data is of dimensions n samplesX m sites
            data = meth_data.data.transpose() # data to test

            # initialize lmm with kinship
            module = lmm.LMM(kinship)
            logging.info('Running LMM...')
            
            t0 = time.time()
            if not args.oneld: # run lmm for each site so logdelta will be calculated for each site (TODO sometime move this option as an argument of LMM class and not of the parser (now, parser calls LMM class with different site each time thats patchy)
              logging.info("LMM calculating logdelta for each site, this will take a while...")
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
                res = module.run(data_site_i, pheno, covars, [meth_data.cpgnames[i]], args.norm, args.reml)
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
              logging.info("computing log delta...")
              lmm_results = module.run(data, pheno, covars, meth_data.cpgnames, args.norm, args.reml)
              cpgnames, pvalues, intercepts_betas, covars_betas, sites_betas, sigmas_e, sigmas_g, stats =  lmm_results
            

            logging.debug("LMM is done in %0.2f seconds" %(time.time()-t0))

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
