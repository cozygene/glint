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
        output_filename = output_file_prefix + self.SCATTER_OUTPUT_FILE
        pca_scatter_plot = plot.PCAScatterPlot(pca_out, plots_number = amount, save_file = output_filename)
        pca_scatter_plot.draw()
        logging.info("pca outliers scatter plot is saved to %s"% output_filename)


    def exclude_maxpcstds(self, pcstds):
        """
        pcstds is a list of lists (or tuples) where the first index is the pc_index and the second index is the std_num
        exclude samples that have std above std_num times higer or lower than the pc std on every std_index
        
        for example, let pcstds be [(1,3),(5,4)] - that will exclude samples that have std above 3 or below -3 on pc 1 and above 4 or below -4 in pc 5 
        """
        pca_out = pca.PCA(self.meth_data.data.transpose())

        remove_samples_indices = set()
        for pc_index, std_num in pcstds:
            logging.info("finding samples with std higher than %d stds or lower than %d stds on pc number %d..." % (std_num, -1 * std_num, pc_index))
            pc = pca_out.P[:,pc_index]
            std_pc = std(pc)
            remove_samples_indices.update(where((pc > std_num * std_pc) | (pc < -1 * std_num * std_pc))[0])

        if remove_samples_indices:
            logging.info("excluding samples with max std...")
            self.meth_data.remove_samples_indices(list(remove_samples_indices))
