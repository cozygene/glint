from pandas import Index, unique
from numpy import loadtxt, empty, where, ndarray
import logging
import os

SITES_INFO_FILE = os.path.join(os.path.dirname(__file__), "assets/HumanMethylationSites")

class SitesInfo(object):

    def __init__(self, cpgnames, chromosomes, positions, genes, categories):
        self.cpgnames = cpgnames
        self.chromosomes = chromosomes #could be int (1,2..) or string ('X', 'Y') therfore - string type
        self.positions = positions.astype(int)
        self.genes = genes
        self.categories = categories

class SitesInfoGenerator(SitesInfo):
    """
    SiteInfo object holds methylation site information for a list of sites
    cpgnames_list - the list of sites to extract and hold the information on.
    """

    CPGNAME_INDEX = 0
    CHROMOSOME_INDEX = 1
    POSITION_INDEX = 2
    GENE_INDEX = 3
    ISLAND_INDEX = 4

    def __init__(self, cpgnames_list):
        logging.info("Loading information about methylation sites...")
        data = loadtxt(SITES_INFO_FILE, dtype = str, ndmin = 2, delimiter=',')
        all_cpgnames = data[:, self.CPGNAME_INDEX]
        indices = self._get_cpgnames_indicis(all_cpgnames, cpgnames_list)
        relevant_data = data[indices, :]
        super(SitesInfoGenerator, self).__init__(cpgnames_list,                        \
                                                 relevant_data[:, self.CHROMOSOME_INDEX], \
                                                 relevant_data[:, self.POSITION_INDEX],   \
                                                 relevant_data[:, self.GENE_INDEX],       \
                                                 relevant_data[:, self.ISLAND_INDEX])
        


    def _get_cpgnames_indicis(self, all_cpgnames, cpgnames_list):
        """
        returns the indexes in all_cpgnames of the cpgnames that found in cpgnames_list.
        the indexes are ordered by cpgnames_list:
            if output is this function output so all_cpgnames[output[0]] is equal to  cpgnames[0]
        """
        logging.info("Searching for relevant methylation sites information...")
        
        sorted_indices = cpgnames_list.argsort()  # the indexes of cpgnames_list when values are sorted
        orderd_indices = empty((len(sorted_indices)), dtype = int)# will hold the indexes in the order of cpgnames_list

        # get the indexes of all_cpgnames where  cpgnames_list found by the order of cpgnames_list 
        assert type(cpgnames_list) == list or type(cpgnames_list) == ndarray
        assert type(all_cpgnames) == list or type(all_cpgnames) == ndarray
        indexes_list = where(Index(unique(cpgnames_list)).get_indexer(all_cpgnames) >= 0)[0] # faster than np.where(np.in1d
        orderd_indices[sorted_indices] = indexes_list
        return orderd_indices