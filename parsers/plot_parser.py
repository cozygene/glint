import os
import sys
import argparse
import logging
from utils import common, plot, pca
from module_parser import ModuleParser
from modules import ewas, kit
from parsers import MethylationDataParser

QQ_PLOT_SUFFIX = ".glint.qqplot"
MANHATTEN_SUFFIX = ".glint.manhattan"

class QQPlotParser(ModuleParser):
    def __init__(self, parser):
        """
        QQPlot Notes:
        can be run with --ewas test or get an ewas test result file (with --result flag):
          example: glint.py --ewas --lmm --plot --qqplot --datafile...
                or glint.py --plot --qqplot --result <file with the output of an ewas test>

        * validates that either a result file was supplied or it was executed with ewas test

        * by default plot will have no title unless user supplies a title with --title flag

        will save output file in both .png and .eps formats
        """
        plot = parser.add_argument_group('qqplot', 'Generates qq-plot')
        plot.add_argument('--results', type = argparse.FileType('r'),  help = "an EWAS test results file. Supply this if --ewas test was not selected") 
        plot.add_argument('--title', type = str,  help = "The title for the plot, will be left empty if not supplied")
        
        super(QQPlotParser, self).__init__(plot)


    def run(self, args, ewas_result_obj = None):
        # not result file nor ewas result supplied
        if ewas_result_obj is None:
            common.terminate("Must supply results file to qq-plot. Use --result to supply a glint results file or --ewas to run a new test")
        
        # plot the p-value
        output_perfix = args.out
        qqplot_out = "results" + QQ_PLOT_SUFFIX if output_perfix is None else output_perfix + QQ_PLOT_SUFFIX
        qqplot = plot.QQPlot(save_file = qqplot_out)
        
        qqplot.draw(ewas_result_obj.pvalues, title = args.title) # todo add option for the user to choose x and y titles


class ManhattanPlotParser(ModuleParser):
    def __init__(self, parser):
        """
        ManhattanPlot Notes:
        can be run with --ewas test or get an ewas test result file (with --result flag):
          example: glint.py --ewas --lmm --plot --manhattan --datafile...
                or glint.py --plot --manhattan --result <file with the output of an ewas test>

        * validates that either a result file was supplied or it was executed with ewas test
        
        * by default plot will have no title unless user supplies a title with --title flag
    
        will save output file in both .png and .eps formats
        """
        plot = parser.add_argument_group('manhattan', 'Plots Manhattan plot')
        plot.add_argument('--results', type = argparse.FileType('r'),  help = "An EWAS test results file. Supply this if --ewas test was not selected")
        plot.add_argument('--title', type = str,  help = "The title for the plot")

        super(ManhattanPlotParser, self).__init__(plot)


    def run(self, args, ewas_result_obj = None):
        # not result file nor ewas result supplied
        if ewas_result_obj is None:
            common.terminate("Must supply results file. Use --result to supply a glint result file or --ewas to run a new test")
        
        # plot the p-value
        output_perfix = args.out
        manplot_out = "results" + MANHATTEN_SUFFIX if output_perfix is None else output_perfix + MANHATTEN_SUFFIX
        manplot = plot.ManhattanPlot(save_file = manplot_out)

        manplot.draw(ewas_result_obj.cpgnames, ewas_result_obj.pvalues, ewas_result_obj.sites_info.chromosomes, ewas_result_obj.sites_info.positions, title = args.title)


class PCAScatterPlotParser(ModuleParser):
  SCATTER_OUTPUT_FILE = "plotpcs"
  def __init__(self, parser):
    """
    PCAScatterPlot Notes:

    runs PCA and plots couples of PCs:
      number of PCs couples to plot specified with --numpcs  flag
      if --numpcs is i there will be (i-1) plots with the following couples: PC(1) vs PC(2), PC(2) vs PC(3),... PC(i-1) vs PC(i)
      * assumes that (numpcs + 1) < (number of samples in data)

    - plot will have no title

    to run:
          glint.py --plot --plotpcs --numpcs 3 --datafile ...


    output file saved in both .png and .eps formats
    """
    pca = parser.add_argument_group('plotpcs', 'Generates scatter plots of the first several PCs of the data')
    pca.add_argument('--numpcs', required = True, type = int, default = 2, help = "Number of PCs to plot") 


    self.meth_data_parser = MethylationDataParser(parser)
    super(PCAScatterPlotParser, self).__init__(pca)

  def validate_args(self, args):
    # add the methylation data parser since user can choose that module flags when he wants to plot PCA scatter
    self.meth_data_parser.validate_args(args)
    self.all_args.extend(self.meth_data_parser.all_args)
    self.required_args.extend(self.meth_data_parser.required_args)

  def run(self, args, meth_data):
    # run pca and plot PCs
    output_filename = args.out if args.out else self.SCATTER_OUTPUT_FILE

    try:
      assert args.numpcs + 1 < meth_data.samples_size

      logging.info("Running PCA...")
      pca_out = pca.PCA(meth_data.data.transpose()) # meth_data should be transposed before passing to pca
      
      logging.info("Plotting first %s PCs..." % args.numpcs)
      pca_scatter_plot = plot.PCAScatterPlot(pca_out, plots_number = args.numpcs, save_file = output_filename)
      pca_scatter_plot.draw()

    except Exception:
      logging.exception("In pca plot parser")
      raise
    


class PlotParser(ModuleParser):
  """
  Plot Notes:
  allows to run only one plot at a time
  plots options:
    - qq plot
    - manhattan plot
    - pca scatter plot

  to see each plot flags and documentation refer to it's parser.

  all plots saves output file in .png and .eps formats.
  """
  def __init__(self, parser):
    plot = parser.add_argument_group('plot', 'For data visualization.')
    plot.add_argument('--qqplot', action='store_true',   help = "QQ-plot")
    plot.add_argument('--manhattan', action='store_true', help = "Manhattan plot")
    plot.add_argument('--plotpcs', action='store_true', help = "PCA scatter plot")
    
    self.qqplot_parser = QQPlotParser(parser)
    self.manhattan_parser = ManhattanPlotParser(parser)
    self.plotpcs_parser = PCAScatterPlotParser(parser)
    super(PlotParser, self).__init__(plot)


  def validate_args(self, args):
    plot_counter = 0

    if args.qqplot:
        plot_counter += 1
        self.qqplot_parser.validate_args(args)
        self.all_args.extend(self.qqplot_parser.all_args)
        self.required_args.extend(self.qqplot_parser.required_args)

    if args.plotpcs:
        plot_counter += 1
        self.plotpcs_parser.validate_args(args)
        self.all_args.extend(self.plotpcs_parser.all_args)
        self.required_args.extend(self.plotpcs_parser.required_args)

    if args.manhattan:
        plot_counter += 1
        self.manhattan_parser.validate_args(args)
        self.all_args.extend(self.manhattan_parser.all_args)
        self.required_args.extend(self.manhattan_parser.required_args)

    if plot_counter == 0:
        common.terminate("Please select a plot type")
    # if plot_counter >1:
    #     common.terminate("plese select only one plot option")

    super(PlotParser, self).validate_args(args)

  def run(self, args, meth_data = None, ewas_result_obj = None):
    """
     assumes or --result is supplied or ewas_result_obj is supplied (not None)
     (this is backed-up in the validate_args) function
    """
    try:
      # both result file and results from ewas test supplied
      if args.results and ewas_result_obj is not None: # user ran ewas and plot together but supplied results file - do not plot cause it could be 
          common.terminate("Couldn't choose between ewas results file %s and the new test results" % args.results.filename) #todo filename /file?

      #supplied result file - extract results
      if args.results:
          ewas_result_obj = ewas.EWASResultsParser(args.results)

      if args.qqplot:
          self.qqplot_parser.run(args, ewas_result_obj)

      if args.plotpcs:
          assert(meth_data)
          self.plotpcs_parser.run(args, meth_data)

      if args.manhattan:    
          self.manhattan_parser.run(args, ewas_result_obj)

    except Exception:
      logging.exception("in plotting")
      raise
    
