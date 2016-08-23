import os
import sys
import argparse
import logging
from utils import common, plot
from module_parser import ModuleParser
from modules import ewas

QQ_PLOT_SUFFIX = ".glint.qqplot"
MANHATTEN_SUFFIX = ".glint.manhattan"

class PlotParser(ModuleParser):

  def __init__(self, parser):
    plot = parser.add_argument_group('plot', 'Plotting options TODO Elior, add description which will be shown when --help')
    plot.add_argument('--result', type = argparse.FileType('r'),  help = "an EWAS test results file (glint format). Supply this if --ewas test was not selected")
    plot.add_argument('--title', type = str,  help = "the title for the plot")
    plot.add_argument('--qqplot', action='store_true',   help = "QQ-plot")
    plot.add_argument('--manhattan', action='store_true', help = "Manhattan plot")
       
    super(PlotParser, self).__init__(plot)


  def validate_args(self, args):
    # if user didnt choose --ewas it means that he won't run EWAS test so he must supply his own results.
    # in such case --result flag is required (so the user can supply the results he wants to plot)
    if not args.ewas:
        self.required_args.append('result')

    plot_counter = 0
    if args.manhattan:
      plot_counter += 1
    if args.qqplot:
      plot_counter +=1

    if plot_counter == 0:
        common.terminate("plese select plot type --qqplot, --manhattan e.g")
    if plot_counter >1:
        common.terminate("plese select only one plot option")
    #todo some default  plot (qq or manhattan)

    super(PlotParser, self).validate_args(args)

  def run(self, args, ewas_result_obj = None):
    """
     assumes or --result is supplied or ewas_result_obj is supplied (not None)
     (this is backed-up in the validate_args) function
    """
    try:
      if args.result and ewas_result_obj is not None: # user ran ewas and plot together but supplied results file - do not plot cause it could be 
        logging.warning("couldn't chose between ewas results file %s and the new test results" % args.result.filename) #todo filename /file?
        return

      if args.result:
        ewas_result_obj = ewas.EWASResultsParser(args.result)
      
      import pdb
      # pdb.set_trace()

      output_perfix = args.out

      if args.qqplot :
        # plot the p-value
        qqplot_out = "results" + QQ_PLOT_SUFFIX if output_perfix is None else output_perfix + QQ_PLOT_SUFFIX
        qqplot = plot.QQPlot(save_file = qqplot_out)
        
        qqplot.draw(ewas_result_obj.pvalues, title = args.title) # todo add option for the user to choose x and y titles

      if args.manhattan:    
        manplot_out = "results" + MANHATTEN_SUFFIX if output_perfix is None else output_perfix + MANHATTEN_SUFFIX
        manplot = plot.ManhattanPlot(save_file = manplot_out)
        manplot.draw(ewas_result_obj.cpgnames, ewas_result_obj.pvalues, ewas_result_obj.sites_info.chromosomes, ewas_result_obj.sites_info.positions, title = args.title)


    except Exception:
      logging.exception("in plotting")
      raise
    
