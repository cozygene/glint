import os
import sys
import argparse
import logging
from utils import common, plot
from module_parser import ModuleParser
from modules import ewas, kit
from parsers import MethylationDataParser

QQ_PLOT_SUFFIX = ".glint.qqplot"
MANHATTEN_SUFFIX = ".glint.manhattan"

class QQPlotParser(ModuleParser):
    def __init__(self, parser):
        plot = parser.add_argument_group('qqplot', 'Plotting options TODO Elior, add description which will be shown when --help')
        plot.add_argument('--result', type = argparse.FileType('r'),  help = "an EWAS test results file (glint format). Supply this if --ewas test was not selected")
        plot.add_argument('--title', type = str,  help = "the title for the plot")
        
        super(QQPlotParser, self).__init__(plot)


    def run(self, args, ewas_result_obj = None):
        # both result file and results from ewas test supplied
        if args.result and ewas_result_obj is not None: # user ran ewas and plot together but supplied results file - do not plot cause it could be 
            logging.warning("couldn't choose between ewas results file %s and the new test results" % args.result.filename) #todo filename /file?
            return

        # not result file nor ewas result supplied
        if ewas_result_obj is None and args.result is None:
            common.terminate("must supply results to qq-plot. use --result to supply a glint result file or --ewas to run a new test")
        
        #supplied result file - extract results
        if args.result:
            ewas_result_obj = ewas.EWASResultsParser(args.result)

        # plot the p-value
        output_perfix = args.out
        qqplot_out = "results" + QQ_PLOT_SUFFIX if output_perfix is None else output_perfix + QQ_PLOT_SUFFIX
        qqplot = plot.QQPlot(save_file = qqplot_out)
        
        qqplot.draw(ewas_result_obj.pvalues, title = args.title) # todo add option for the user to choose x and y titles


class ManhattanPlotParser(ModuleParser):
    def __init__(self, parser):
        plot = parser.add_argument_group('manhattan', 'Plotting options TODO Elior, add description which will be shown when --help')
        plot.add_argument('--result', type = argparse.FileType('r'),  help = "an EWAS test results file (glint format). Supply this if --ewas test was not selected")
        plot.add_argument('--title', type = str,  help = "the title for the plot")

        super(ManhattanPlotParser, self).__init__(plot)


    def run(self, args, ewas_result_obj = None):
        # both result file and results from ewas test supplied
        if args.result and ewas_result_obj is not None: # user ran ewas and plot together but supplied results file - do not plot cause it could be 
            logging.warning("couldn't choose between ewas results file %s and the new test results" % args.result.filename) #todo filename /file?
            return

        # not result file nor ewas result supplied
        if ewas_result_obj is None and args.result is None:
            common.terminate("must supply results to manhattan-plot. use --result to supply a glint result file or --ewas to run a new test")
        
        #supplied result file - extract results
        if args.result:
            ewas_result_obj = ewas.EWASResultsParser(args.result)

        # plot the p-value
        output_perfix = args.out
        manplot_out = "results" + MANHATTEN_SUFFIX if output_perfix is None else output_perfix + MANHATTEN_SUFFIX
        manplot = plot.ManhattanPlot(save_file = manplot_out)

        manplot.draw(ewas_result_obj.cpgnames, ewas_result_obj.pvalues, ewas_result_obj.sites_info.chromosomes, ewas_result_obj.sites_info.positions, title = args.title)


class PCAScatterPlotParser(ModuleParser):#MethylationDataParser):
  """
  All the tools provided to the user which are not a module (don't have a run() function)
  """
  def __init__(self, parser):
    pca = parser.add_argument_group('plotpcs', 'TODO Elior, add description which will be shown wjen --help')
    pca.add_argument('--numpcs', required = True, type = int, help = " number of pcs to plot TODO Elior, edit")

    self.meth_data_parser = MethylationDataParser(parser)
    super(PCAScatterPlotParser, self).__init__(pca)

  def validate_args(self, args):
    # if user didnt choose --ewas it means that he won't run EWAS test so he must supply his own results.
    # in such case --result flag is required (so the user can supply the results he wants to plot)
    plot_counter = 0

    self.meth_data_parser.validate_args(args)
    self.all_args.extend(self.meth_data_parser.all_args)

  def run(self, args, meth_data):
    output_perfix = args.out

    try:
      self.pca_utils = kit.PCAKit(meth_data)

      assert args.numpcs + 1 < meth_data.samples_size
      self.pca_utils.draw_pca_scatter(args.numpcs, output_perfix)

    except Exception:
      logging.exception("in utils")
      raise
    


class PlotParser(ModuleParser):

  def __init__(self, parser):
    plot = parser.add_argument_group('plot', 'Plotting options TODO Elior, add description which will be shown when --help')
    plot.add_argument('--qqplot', action='store_true',   help = "QQ-plot")
    plot.add_argument('--manhattan', action='store_true', help = "Manhattan plot")
    plot.add_argument('--plotpcs', action='store_true', help = "PCA scatter plot")
    
    self.qqplot_parser = QQPlotParser(parser)
    self.manhattan_parser = ManhattanPlotParser(parser)
    self.plotpcs_parser = PCAScatterPlotParser(parser)
    super(PlotParser, self).__init__(plot)


  def validate_args(self, args):
    # if user didnt choose --ewas it means that he won't run EWAS test so he must supply his own results.
    # in such case --result flag is required (so the user can supply the results he wants to plot)
    plot_counter = 0

    if args.qqplot:
        plot_counter += 1
        self.qqplot_parser.validate_args(args)
        self.all_args.extend(self.qqplot_parser.all_args)

    if args.plotpcs:
        plot_counter += 1
        self.plotpcs_parser.validate_args(args)
        self.all_args.extend(self.plotpcs_parser.all_args)

    if args.manhattan:
        plot_counter += 1
        self.manhattan_parser.validate_args(args)
        self.all_args.extend(self.manhattan_parser.all_args)

    if plot_counter == 0:
        common.terminate("plese select plot type")
    if plot_counter >1:
        common.terminate("plese select only one plot option")

    super(PlotParser, self).validate_args(args)

  def run(self, args, meth_data = None, ewas_result_obj = None):
    """
     assumes or --result is supplied or ewas_result_obj is supplied (not None)
     (this is backed-up in the validate_args) function
    """
    try:
      # if args.result and ewas_result_obj is not None: # user ran ewas and plot together but supplied results file - do not plot cause it could be 
      #   logging.warning("couldn't choose between ewas results file %s and the new test results" % args.result.filename) #todo filename /file?
      #   return

      # if args.result:
      #   ewas_result_obj = ewas.EWASResultsParser(args.result)
      
      # output_perfix = args.out

      if args.qqplot :
          return self.qqplot_parser.run(args, ewas_result_obj)
      #   # plot the p-value
      #   qqplot_out = "results" + QQ_PLOT_SUFFIX if output_perfix is None else output_perfix + QQ_PLOT_SUFFIX
      #   qqplot = plot.QQPlot(save_file = qqplot_out)
        
      #   qqplot.draw(ewas_result_obj.pvalues, title = args.title) # todo add option for the user to choose x and y titles

      if args.plotpcs:
          assert(meth_data)
          return self.plotpcs_parser.run(args, meth_data)

      if args.manhattan:    
          return self.manhattan_parser.run(args, ewas_result_obj)
      #   manplot_out = "results" + MANHATTEN_SUFFIX if output_perfix is None else output_perfix + MANHATTEN_SUFFIX
      #   manplot = plot.ManhattanPlot(save_file = manplot_out)
      #   manplot.draw(ewas_result_obj.cpgnames, ewas_result_obj.pvalues, ewas_result_obj.sites_info.chromosomes, ewas_result_obj.sites_info.positions, title = args.title)


    except Exception:
      logging.exception("in plotting")
      raise
    
