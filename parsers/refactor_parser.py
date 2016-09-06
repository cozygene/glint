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
from numpy import loadtxt

# list of files containing bad sites ids
BAD_PROBES_FILES = [
                    "assets/48639-non-specific-probes-Illumina450k.txt",
                    "assets/artifacts_chen.2013.txt",
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
      --minstd: minstd cannot be greater than 1 and smaller than 0. program terminates otherwise.
      --numcomp: number of components must be at least k and smaller than the number of samples size. program terminates otherwise.


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
      refactor.add_argument('--fs',      type = str, default = 'normal', help = "feature selection mode; options: normal, controls, phenotype (DEFAULT=normal)")
      refactor.add_argument('--minstd',  type = float, default = 0.02, help = "threshold for excluding low variance sites (DEFAULT=0.02) (all sites with std lower than this threshold will be excluded)") # all sites with std lower than --minstd will be excluded 
      refactor.add_argument('--rmcovar', action = "store_true", help = "specify if you want to remove covariates from data") # if called with --rmcovar flag, the covariates are  removed from the data. (by default (without the flag) they are not removed). Note that this flag appears also in epistructure

      super(RefactorParser, self).__init__(refactor)


    def run(self, args, meth_data, output_perfix = None):
      try:
        if output_perfix is None:
          output_perfix = "output"
        bad_probes_list = set()
        [bad_probes_list.update(loadtxt(os.path.join( os.path.dirname(__file__), probes_file), dtype=str)) for probes_file in BAD_PROBES_FILES]
        self.module  = refactor.Refactor(methylation_data = meth_data, 
                              k = args.k, 
                              t = args.t, 
                              minstd = args.minstd,
                              feature_selection = args.fs.lower().strip(), 
                              num_components = args.numcomp,
                              remove_covars = args.rmcovar,
                              bad_probes_list = bad_probes_list,
                              ranked_output_filename = output_perfix + "." + refactor.RANKED_FILENAME, 
                              components_output_filename  = output_perfix + "." + refactor.COMPONENTS_FILENAME)
        self.module.run()
      except Exception :
        logging.exception("in refactor")
        raise
      
