import sys
import logging
from utils import common
import argparse
# from module_parser import ModuleParser
# from modules import predictor
from pickle import load

                
# snps_info = load(open("parsers/assets/model-snps_info.pickle", 'r'))
sites_info = load(open("parsers/assets/model-sites_info.pickle", 'r'))

import pdb
pdb.set_trace()

# def get_model_without_artifacts():
#     blaa
# class PredictorParser(ModuleParser):
#     def __init__(self, parser):

#         predictor = parser.add_argument_group('predictor', 'TODO Elior,add ewas description here')
#         predictor.add_argument('--plink',  type = argparse.FileType('r'), required = True, help = "predict methylation level by SNPs information in plink file")
#         predictor.add_argument('--score',  type = float, default = 0.5, help = "a score S(between 0 and 1) - predict only cpgs with prediction correlation score at least S ")
#         super(PredictorParser, self).__init__(predictor)
        

#     def run(self, args, output_perfix = ""):
#         try:
#             self.module  = predictor.Predictor(get_model_without_artifacts())
#             self.module.run()
#         except Exception :
#             logging.exception("in predictor")
#             raise