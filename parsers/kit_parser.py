import os
import sys
import argparse
import logging
from modules import kit
from utils import common
from module_parser import ModuleParser

class KitParser(ModuleParser):
  """
  All the tools provided to the user which are not a module (don't have a run() function)
  """
  def __init__(self, parser):
    pca = parser.add_argument_group('pca', 'todo')
    pca.add_argument('--maxpcstd', type = int, action = 'append', nargs = 2, help = "todo pc index and std number of times for removing outliers")
    pca.add_argument('--plotpcs', type = int, help = "todo number of pcs to plot")

    super(KitParser, self).__init__(pca)

  def run(self, args, meth_data, output_perfix = ""):
    try:
      self.pca_utils = kit.PCAKit(meth_data)
      if args.maxpcstd is not None:
        self.pca_utils.exclude_maxpcstds(args.maxpcstd)
      if args.plotpcs is not None:
        assert args.plotpcs + 1 < meth_data.samples_size
        self.pca_utils.draw_pca_scatter(args.plotpcs)


    except Exception:
      logging.exception("in utils")
      raise
    
