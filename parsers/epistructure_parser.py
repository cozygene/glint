import sys
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import epistructure
from numpy import loadtxt
                
INFORMATIVE_ANCESTRY_CPG_LIST = "parsers/assets/KORA_model_multiple_snps_W_50_M_10_cpgs_with_corr_th_0.5.txt"
DEAFUALT_FILENAME = 'ancestry_pcs'
class EpistructureParser(ModuleParser):
    def __init__(self, parser):

        epistructure_parser = parser.add_argument_group('epistructure', 'TODO Elior,add epistructure description here')
        epistructure_parser.add_argument('--savepcs', type = int, default = 2, help = "number of captured ancestry PCs to save TODO Elior, edit")
        epistructure_parser.add_argument('--rmcovar', action = "store_true", help = "specify if you want to remove covariates from data") # if epistructure is called with --rmcovar flag, the covariates are  removed from the data. (by default (without the flag) they are not removed) 
        super(EpistructureParser, self).__init__(epistructure_parser)
        

    def run(self, args, meth_data, output_perfix = None):
        if output_perfix is None:
            output_perfix = DEAFUALT_FILENAME
        try:
            informative_sites = loadtxt(INFORMATIVE_ANCESTRY_CPG_LIST, dtype = str)
            self.module = epistructure.Epistructure(meth_data, informative_sites)
            self.module.capture_ancestry(args.savepcs, args.rmcovar, output_perfix)
        except Exception :
            logging.exception("in epistructure")
            raise



