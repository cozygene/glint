"""""
ReFactor
"""""
import os
import sys
import argparse
import logging
from modules import refactor
from utils import common
from module_parser import ModuleParser
from numpy import loadtxt, array

# list of files containing bad sites ids
BAD_PROBES_FILES = [
                    os.path.join(os.path.dirname(__file__),"assets/HumanMethylationSites_X_Y.txt"),
                    os.path.join(os.path.dirname(__file__), "assets/nonspecific_probes.txt"), # those sites might have large variance , so should be removed
                    os.path.join(os.path.dirname(__file__), "assets/polymorphic_cpgs.txt")
                   ]

class RefactorParser( ModuleParser ):

    def __init__(self, parser):
      """
      Flags:
      --fs feature selection options are normal, controls (validates that phenotype must is binary), phenotype (validates that phenotype file was  supplied)
      --k: k must be at least 2 and smaller than the number of samples size.
             if controls feature selected -  k cannot be greater than amount of controls samples found in the controls phenotype.
              program terminates otherwise.
      --t: t cannot be greater than the number of sites or smaller than k. program terminates otherwise.
      --stdth:  cannot be greater than 1 and smaller than 0. program terminates otherwise.
      --numcomp: number of components must be at least k and smaller than the number of samples size. program terminates otherwise.
      --covar - list of the names of covariates to use.
                if flag is not set will not use any covariates. if set and no name if specified will use all of the covariates. o
                therwise (names specified) will only use the covariates which name is in the list
      --pheno - list of names of phenotypes to use. behaves liske --covar

      ReFactor Notes:
       - removes bad sites by default (and there is no option to disable it). the files with the bad sites list are found in BAD_PROBES_FILES 
       - removes sites with low std (lower than the threshold at --minstd)
       - doesn't remove covariates by default (use --rmcovar to remove them)
       - output:
          - saves eanked sites list
          - save the ReFactor components


      """
      refactor = parser.add_argument_group('refactor', 'Additional options if --refactor is selected')

      refactor.add_argument('--k',       type = int, required = True, help = "The number of assumed cell types")
      refactor.add_argument('--t',       type = int, default = 500, help = "The number of sites to use for computing the ReFACtor components (DEFAULT=500)")
      refactor.add_argument('--numcomp', type = int, help = "The number of ReFACTor components to output (DEFAULT=K)")
      refactor.add_argument('--fs',      type = str, default = 'normal', help = "Feature selection mode; options: normal, controls, phenotype (DEFAULT=normal)")
      refactor.add_argument('--stdth',  type = float, default = 0.02, help = "Threshold for excluding low variance sites (DEFAULT=0.02) (all sites with std lower than this threshold will be excluded)") 
      refactor.add_argument('--covar', type = str, nargs='*', help = "List of covariate names to use.")
      refactor.add_argument('--pheno', type = str, nargs='*', help = "The phenotype name to use (if 'phenotype' was selected under --fs)")
        
      super(RefactorParser, self).__init__(refactor)


    def run(self, args, meth_data, output_perfix = None):
      try:
        if args.pheno is not None and meth_data.phenotype is None:
          common.terminate("There is no phenotype in the data, use --phenofile to add phenotype.")
        if not output_perfix:
          output_perfix = "output"
        bad_probes_list = set()
        [bad_probes_list.update(loadtxt(probes_file, dtype=str)) for probes_file in BAD_PROBES_FILES]
        bad_probes_list = array(list(bad_probes_list))
        self.module  = refactor.Refactor(methylation_data = meth_data, 
                              k = args.k, 
                              t = args.t, 
                              stdth = args.stdth,
                              feature_selection = args.fs.lower().strip(), 
                              num_components = args.numcomp,
                              use_covars = args.covar,
                              use_phenos = args.pheno,
                              bad_probes_list = bad_probes_list,
                              ranked_output_filename = output_perfix + "." + refactor.RANKED_FILENAME, 
                              components_output_filename  = output_perfix + "." + refactor.COMPONENTS_FILENAME)
        self.module.run()
      except Exception :
        logging.exception("in refactor")
        raise
      
