import os
import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import epistructure
from numpy import loadtxt

# list of sites that are ancestry-informative                
INFORMATIVE_ANCESTRY_CPG_LIST = os.path.join(os.path.dirname(__file__), "assets/epistructure_reference_sites.txt")

class EpistructureParser(ModuleParser):
    def __init__(self, parser):
        """
        --savepcs:   number of captured ancestry PCs to save
        """
        epistructure_parser = parser.add_argument_group('epistructure', 'The EPISTRUCTURE algorithm captures population structure from methylation data')
        epistructure_parser.add_argument('--savepcs', type = int, default = 1, help = "Number of EPISTRUCTURE PCs to save (DEFAULT=1)")
        epistructure_parser.add_argument('--covar', type = str, nargs='*', help = "List of covariate names to use")
      
        super(EpistructureParser, self).__init__(epistructure_parser)
        

    def run(self, args, meth_data, output_perfix = None):
        output_filename = epistructure.EPISTRUCTURE_FILE_SUFFIX if output_perfix is None else output_perfix + "." +  epistructure.EPISTRUCTURE_FILE_SUFFIX
        try:
            informative_sites = loadtxt(INFORMATIVE_ANCESTRY_CPG_LIST, dtype = str)
            informative_sites = common.loadtxt(INFORMATIVE_ANCESTRY_CPG_LIST, dtype = str)
            self.module = epistructure.Epistructure(meth_data, informative_sites)
            self.module.capture_ancestry(args.savepcs, args.covar, output_filename)
            return self.module.components
        except Exception :
            logging.exception("in epistructure")
            raise



