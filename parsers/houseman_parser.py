import sys
import os
import logging
from utils import common
import argparse
from module_parser import ModuleParser
from modules import houseman

HOUSEMAN_DEFAULT_REFERENCE = os.path.join(os.path.dirname(__file__), "assets/12859_2016_943_MOESM5_ESM_ref.txt")
HOUSEMAN_OUTPUT_NAME = "houseman_estimates.txt"

class HousemanParser(ModuleParser):
    def __init__(self, parser):
        houseman = parser.add_argument_group('houseman todo elior add description') # numbering in the group name because help print it by abc order
        houseman.add_argument('--reference', type = argparse.FileType('r'), default = open(HOUSEMAN_DEFAULT_REFERENCE, 'r'),  help = "reference file in a specific format. see documentation for more details")
        super(HousemanParser, self).__init__(houseman)

    def run(self, args, meth_data, output_perfix):
        logging.info("running Houseman...")
        filename = HOUSEMAN_OUTPUT_NAME
        if output_perfix:
            filename = output_perfix + "." + filename
        
        self.module = houseman.Houseman(meth_data, args.reference, filename)
