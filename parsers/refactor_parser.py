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

BAD_PROBES_FILES = [
                    "assets/48639-non-specific-probes-Illumina450k.txt",
                    "assets/artifacts_chen.2013.txt",
                   ]

class RefactorParser( ModuleParser ):

    def __init__(self, parser):
      refactor = parser.add_argument_group('refactor', 'Additional options if --refactor is selected')

      refactor.add_argument('--k',       type = int, required = True, help = "The number of assumed cell types")
      refactor.add_argument('--t',       type = int, default = 500, help = "The number of sites to use for computing the ReFACtor components (DEFAULT=500)")
      refactor.add_argument('--numcomp', type = int, help = "The number of ReFACTor components to output (DEFAULT=K)")
      refactor.add_argument('--fs',      type = str, default = 'normal', help = "feature selection mode; options: normal, controls, phenotype (DEFAULT=normal)")
      refactor.add_argument('--minstd',   type = float, default = 0.02, help = "threshold for excluding low variance sites")
      refactor.add_argument('--suppress_covars', type = bool, default = False, help = "set to True if you don't need to remove covariates ")

      super(RefactorParser, self).__init__(refactor)


    def run(self, args, meth_data, output_perfix = ""):
      try:

        bad_probes_list = set()
        [bad_probes_list.update(loadtxt(probes_file, dtype=str)) for probes_file in BAD_PROBES_FILES]
        self.module  = refactor.Refactor(methylation_data = meth_data, 
                              k = args.k, 
                              t = args.t, 
                              minstd = args.minstd,
                              feature_selection = args.fs.lower().strip(), 
                              num_components = args.numcomp,
                              suppress_covars = args.suppress_covars,
                              bad_probes_list = bad_probes_list,
                              ranked_output_filename = output_perfix + refactor.RANKED_FILENAME, 
                              components_output_filename  = output_perfix + refactor.COMPONENTS_FILENAME)
        self.module.run()
      except Exception :
        logging.exception("in refactor")
        raise
      
