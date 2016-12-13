from pandas import Index, unique
from numpy import loadtxt, empty, where, ndarray
import logging
import os

SITES_INFO_FILE = os.path.join(os.path.dirname(__file__), "assets/HumanMethylationSites")

class SitesInfo(object):

    def __init__(self, cpgnames, chromosomes, positions, genes, categories):
        self.cpgnames = cpgnames.astype(str)
        self.chromosomes = chromosomes.astype(str) #could be int (1,2..) or string ('X', 'Y') therfore - string type
        self.positions = positions.astype(int)
        self.genes = genes.astype(str)
        self.categories = categories.astype(str)

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
        cpg_found_index, cpg_ref_index = self._get_cpgnames_indicis(all_cpgnames, cpgnames_list)
        
        chromosomes = empty((len(cpgnames_list)), dtype = data.dtype)
        chromosomes.fill("")
        chromosomes[cpg_found_index] = data[cpg_ref_index, self.CHROMOSOME_INDEX]
        
        positions = empty((len(cpgnames_list)), dtype = data.dtype)
        positions.fill("")
        positions[cpg_found_index] = data[cpg_ref_index, self.POSITION_INDEX]
        
        genes = empty((len(cpgnames_list)), dtype = data.dtype)
        genes.fill("")
        genes[cpg_found_index] = data[cpg_ref_index, self.GENE_INDEX]

        categories = empty((len(cpgnames_list)), dtype = data.dtype)
        categories.fill("")
        categories[cpg_found_index] = data[cpg_ref_index, self.ISLAND_INDEX]

        super(SitesInfoGenerator, self).__init__(cpgnames_list,                        \
                                                 chromosomes, \
                                                 positions,   \
                                                 genes,       \
                                                 categories)
        


    def _get_cpgnames_indicis(self, all_cpgnames, cpgnames_list):
        """
        returns the indexes in all_cpgnames of the cpgnames that found in cpgnames_list.
        the indexes are ordered by cpgnames_list:
            if output is this function output so all_cpgnames[output[0]] is equal to  cpgnames[0]
        """
        logging.info("Searching for relevant methylation sites information...")
        
        # get the indexes of all_cpgnames where  cpgnames_list found by the order of cpgnames_list 
        assert type(cpgnames_list) == list or type(cpgnames_list) == ndarray
        assert type(all_cpgnames) == list or type(all_cpgnames) == ndarray
        
        indexes_list = Index(unique(all_cpgnames)).get_indexer(cpgnames_list)  # faster than np.where(np.in1d
        where_found = where(indexes_list>=0)[0]
        what_found = indexes_list[where_found]

        return where_found, what_found
