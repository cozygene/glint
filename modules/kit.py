import logging
from numpy import ceil, sqrt, std, where
from utils import pca, plot
from module import Module



class PCAKit(Module):
    SCATTER_OUTPUT_FILE = "pca_std_outliers_scatter"
    def __init__(self, methylation_data):
        self.meth_data = methylation_data

    def draw_pca_scatter(self, amount, output_file_prefix=''):
        """
        assumes  amount + 1 < self.meth_data.samples_size
        """
        logging.info("running PCA...")
        pca_out = pca.PCA(self.meth_data.data.transpose()) # meth_data should be transposed before passing to pca
        if output_file_prefix is None:
            output_file_prefix = ''
        output_filename = output_file_prefix + self.SCATTER_OUTPUT_FILE
        pca_scatter_plot = plot.PCAScatterPlot(pca_out, plots_number = amount, save_file = output_filename)
        pca_scatter_plot.draw()
