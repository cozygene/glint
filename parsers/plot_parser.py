import os
import sys
import argparse
import logging
from module_parser import ModuleParser

class PlotParser(ModuleParser):

  def __init__(self, parser):
    pca = parser.add_argument_group('plot', 'Plotting options TODO Elior, add description which will be shown wjen --help')
    pca.add_argument('--results', type = argparse.FileType('r'),  help = "A file with glint EWAS test results. TODO Elior, edit number of pcs to plot")
    ewas.add_argument('--qqplot', action='store_const', const='linear_regression',   help = "Run a linear regression analysis (executed by default if --ewas is selected)")
    ewas.add_argument('--manhattan', action='store_const', const='logistic_regression', help = "Run a logistic regression analysis")
       
    super(PlotParser, self).__init__(pca)


  def validate_args(self, args):
    # if user didnt choose --ewas it means that he won't run EWAS test so he must supply his own results.
    # in such case --result flag is required (so the user can supply the results he wants to plot)
    if not args.ewas:
        self.required_args.append('results')

    #todo some default  plot (qq or manhattan)

    super(LMMParser, self).validate_args(args)

  def run(self, args, results_matrix = None):
    """
     assumes or --results is supplied or results_matrix is supplied (not None)
     (this is backed-up in the validate_args) function
    """
    # if results_matrix is None
    # if args.results TODO
    
    except Exception:
      logging.exception("in plotting")
      raise
    
