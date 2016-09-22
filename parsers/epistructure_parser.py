import os
import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import epistructure
from numpy import loadtxt
                
INFORMATIVE_ANCESTRY_CPG_LIST = os.path.join(os.path.dirname(__file__), "assets/epistructure_reference_sites.txt")

class EpistructureParser(ModuleParser):
    def __init__(self, parser):
        """
        --savepcs:   number of captured ancestry PCs to save
        """
        epistructure_parser = parser.add_argument_group('epistructure', 'TODO Elior,add epistructure description here')
        epistructure_parser.add_argument('--savepcs', type = int, default = 1, help = "number of captured ancestry PCs to save TODO Elior, edit")
        epistructure_parser.add_argument('--covar', type = str, nargs='*', help = "list of covariates names to use. If no name is specified will use all the covariates. If flag is not set, will not use any covariate")
      
        super(EpistructureParser, self).__init__(epistructure_parser)
        

    def run(self, args, meth_data, output_perfix = None):
        output_filename = epistructure.EPISTRUCTURE_FILE_SUFFIX if output_perfix is None else output_perfix + "." +  epistructure.EPISTRUCTURE_FILE_SUFFIX
        try:
            informative_sites = loadtxt(INFORMATIVE_ANCESTRY_CPG_LIST, dtype = str)
            self.module = epistructure.Epistructure(meth_data, informative_sites)
            self.module.capture_ancestry(args.savepcs, args.covar, output_filename)
            return self.module.components
        except Exception :
            logging.exception("in epistructure")
            raise



