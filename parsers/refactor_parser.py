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

class RefactorParser( ModuleParser ):

    def __init__(self, parser):
      refactor = parser.add_argument_group('Additional options if --refactor is selected','\n')

      refactor.add_argument('--k',       type = int, required = True, help = "The number of assumed cell types")
      refactor.add_argument('--t',       type = int, default = 500, help = "The number of sites to use for computing the ReFACtor components (DEFAULT=500)")
      refactor.add_argument('--numcomp', type = int, help = "The number of ReFACTor components to output (DEFAULT=K)")
      refactor.add_argument('--fs',      type = str, default = 'normal', help = "feature selection mode; options: normal, controls, phenotype (DEFAULT=normal)")
      refactor.add_argument('--covar',   type = argparse.FileType('r'), help = "A covariates file")
      refactor.add_argument('--pheno',   type = argparse.FileType('r'), help = "A phenotype file")

      super(RefactorParser, self).__init__(refactor)


    def run(self, args, meth_data, output_perfix = ""):
      try:
        self.module  = refactor.Refactor(methylation_data = meth_data, 
                              k = args.k, 
                              t = args.t, 
                              feature_selection = args.fs.lower().strip(), 
                              num_components = args.numcomp, 
                              phenofile = args.pheno,
                              covar = args.covar,
                              ranked_output_filename = output_perfix + refactor.RANKED_FILENAME, 
                              components_output_filename  = output_perfix + refactor.COMPONENTS_FILENAME)
        self.module.run()
      except Exception :
        logging.exception("in refactor")
        raise
      
