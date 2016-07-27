import os
import sys
import copy
import logging
from pickle import dump
from numpy import delete, isnan, nanstd, where, column_stack, std, array, savetxt, in1d
from numpy.ma import average, masked_array
from module import Module
from utils import common, pca
from bisect import bisect_right

COMPRESSED_FILENAME = "methylation_data"
GLINT_FORMATTED_EXTENSION = ".glint" #TODO move to a config file

def validate_no_missing_values(data):
    """
    nan are not supported for version 1.0
    Note: data must be number and not string
    """
    if  isnan(data).sum() > 0:
        common.terminate("missing values are not supported at this version")


class MethylationData(Module):
    """
    meth_data is an numpy matrix
    samples_list, sites_list are numpy arrays
    """
    def __init__(self, meth_data, samples_list, sites_list, phenotype = None ,covar = None):
        self.data = meth_data
        self.samples_ids = samples_list
        self.cpgnames = sites_list
        self.sites_size, self.samples_size = self.data.shape

        if (len(self.samples_ids) != self.samples_size):
            common.terminate("got data with %s samples but %s samples ids" % (self.samples_size, len(self.samples_ids)))


        if (len(self.cpgnames) != self.sites_size):
            common.terminate("got data with %s sites but %s cpgnames" % (self.sites_size, len(self.cpgnames)))

        logging.debug("got methylation data with %s sites and %s samples id" % (self.sites_size, self.samples_size))

        self.phenotype = phenotype
        self.covar = covar


    def exclude_sites_indices(self, sites_indicies_list):
        """
        sites_indicies_list - a list of sites indices to remove from data
        this function removes from the data the cpg sites which indices found in sites_indicies_list
        it updates the sites_size, the cpgnames list
        """
        if len(sites_indicies_list) == 0:
            logging.warning("found no sites to remove")
        else:
            if len(set(range(self.sites_size)).difference(set(sites_indicies_list))) == 0:
                logging.error("all sites are about to be remove") # TODO Elior, should we terminate or warn or error witout terminateing?

            self.data = delete(self.data, sites_indicies_list, axis = 0)
            self.cpgnames = delete(self.cpgnames, sites_indicies_list)
            size_before = self.sites_size
            self.sites_size = len(self.cpgnames)
            logging.debug("%s sites out of %s were excluded (left %s sites)" % (len(sites_indicies_list), size_before, self.sites_size))

    def remove_samples_indices(self, indices_list):
        """
        indices_list - a list of sample indices to remove from data
        removes from the data the samples which indices found in indices_list
        it updates the samples_size, the sample_ids lists
        it also remove the sample from phenotype and covariates file if provided
        """
        if len(indices_list) == 0:
            logging.warning("found no samples to remove")
        else:
            if len(set(range(self.samples_size)).difference(set(indices_list))) == 0:
                logging.error("all samples are about to be remove") # TODO Elior, should we terminate or warn or error witout terminateing?

            self.data = delete(self.data, indices_list, axis = 1)
            if self.phenotype is not None:
                self.phenotype = delete(self.phenotype, indices_list, axis = 0)
            if self.covar is not None:
                self.covar = delete(self.covar, indices_list, axis = 0)
            self.samples_ids = delete(self.samples_ids, indices_list)
            size_before = self.samples_size
            self.samples_size = len(self.samples_ids)
            logging.debug("%s samples out of %s were removed (left %s samples)" % (len(indices_list), size_before, self.samples_size))

    def get_mean_per_site(self):
        """
        returns array that contains the mean value (average) for each methylation site
        """
        masked_data = masked_array(self.data, isnan(self.data)) # create masked array
        return average(masked_data, axis=1)

    def include(self, include_list):
        """
        include_list - a list with sites name to keep in data (remove all other sites)
        this function removes the cpg sites not found in include_list list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("including sites...")
        indices_list = where(in1d(self.cpgnames , include_list))[0]
        logging.info("include sites: %s CpGs from the reference list of %s CpGs were used" % (len(indices_list), len(include_list)))
        self.data = self.data[indices_list, :]
        self.cpgnames = self.cpgnames[indices_list]
        size_before = self.sites_size
        self.sites_size = len(self.cpgnames)
        if (self.data.shape[0] != self.sites_size):
            common.terminate("After including sites, methylation data sites size is %s but we got %s" % (self.data.shape[0], self.sites_size))
        logging.debug("%s sites out of %s were included" % (self.sites_size, size_before))
        logging.info("methylation data new size is %s sites by %s samples" % self.data.shape)


    def exclude(self, exclude_list):
        """
        exclude_list - a list with sites name to remove from data
        this function removes the cpg sites found in self.exclude list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("excluding sites...")
        indices_list = where(in1d(self.cpgnames , exclude_list))[0]
        logging.info("exclude sites: %s CpGs from the reference list of %s CpGs were used" % (len(indices_list), len(exclude_list)))
        self.exclude_sites_indices(indices_list)
        logging.debug("methylation data new size is %s sites by %s samples" % self.data.shape)

    def keep(self, keep_list):
        """
        keep_list- a list of samples ids to keep in data (remove all other sample ids)
        this function removes the samples ids not found in keep_list list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("keeping samples...")
        remove_indices_list = where(False == in1d(self.samples_ids , keep_list))[0]
        logging.info("keep samples: %s samples from the reference list of %s samples were used" % (len(remove_indices_list), len(keep_list)))
        self.remove_samples_indices(remove_indices_list)
        logging.debug("methylation data new size is %s sites by %s samples" % self.data.shape)

    def remove(self, remove_list):
        """
        remove_list - a list of samples ids to remove from data
        this function removes the samples ids found in remove_list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("removing samples...")
        indices_list = where(in1d(self.samples_ids , remove_list))[0]
        logging.info("remove samples: %s samples from the reference list of %s samples were used" % (len(indices_list), len(remove_list)))
        self.remove_samples_indices(indices_list)
        logging.debug("methylation data new size is %s sites by %s samples" % self.data.shape)
        
    def exclude_sites_with_low_mean(self, min_value):
        """
        removes sites with mean < min_value
        """
        logging.info("excluding sites with mean lower than %s..." % min_value)
        min_values_indices = where(self.get_mean_per_site() < min_value)[0]  
        self.exclude_sites_indices(min_values_indices)

    def exclude_sites_with_high_mean(self, max_value):
        """
        removes sites with mean > max_value
        """
        logging.info("removing sites with mean greater than %s..." % max_value)
        max_values_indices = where(self.get_mean_per_site() > max_value)[0]
        self.exclude_sites_indices(max_values_indices)

    def save(self, prefix = ''):
        """
        serializes this object and saves it to methylation_data_filename
        assumes that methylation_data_filename is a valid file 
        """
        if prefix is None:
            prefix = ''
        if prefix != '' and not prefix.endswith('_'):
            prefix = prefix + "_"
        filename = prefix + COMPRESSED_FILENAME
        methylation_data_filename = filename + GLINT_FORMATTED_EXTENSION

        with open(methylation_data_filename, 'wb') as f:
            logging.info("Saving methylation data as glint format at %s" % methylation_data_filename)
            dump(self, f)
        f.close()
        
        logging.info("Saving cpg names to %s" % filename + "_sites_list.txt")
        savetxt(filename + "_sites_list.txt", self.cpgnames, fmt = '%s25')


        logging.info("Saving samples ids to %s" % filename + "_sampless_list.txt")
        savetxt(filename + "_sampless_list.txt", self.samples_ids, fmt = '%s25')

    def remove_lowest_std_sites(self, lowest_std_th = 0.02):
        """
        input: lowest_std_th threshold for excluding low variance sites, all sites with std lower than lowest_std_th will be excluded 
        lowest_std_th is float between 0 and 1
        """
        logging.info("excluding site with variance lower than %s..." % lowest_std_th)
        # get std for each site
        sites_std = nanstd(self.data, axis=1) # calc variance consider NaN
        # sort std and get sites index for each std (sorted, so indices of the lowest std sites will be to the left) 
        std_sorted_indices = sites_std.argsort() # sort the sites_variance and return an array that holds the indices of the sorted values
        # get std list sorted
        std_sorted = sites_std[std_sorted_indices]
        # get the first index in the sorted list which have std higher than lowest_std_th and include all indices started from it
        include_from_index = bisect_right(std_sorted, lowest_std_th)
        if (include_from_index == self.sites_size):
            common.terminate("the provided stdth parameter excludes all sites (stdth = %s)" % lowest_std_th)
        if (include_from_index == 0):
            logging.warning("the provided stdth parameter excludes no sites (stdth = %s)" % lowest_std_th)
        exclude_sites_indices = std_sorted_indices[:include_from_index]

        # exclude all sites with low std
        self.exclude_sites_indices(std_sorted_indices[:include_from_index])

    def exclude_maxpcstds(self, pcstds):
        """
        pcstds is a list of lists (or tuples) where the first index is the pc_index and the second index is the std_num
        exclude samples that have std above std_num times higer or lower than the pc std on every std_index
        
        for example, let pcstds be [(1,3),(5,4)] - that will exclude samples that have std above 3 or below -3 on pc 1 and above 4 or below -4 in pc 5 
        """
        pca_out = pca.PCA(self.data.transpose())

        maxpcstds_samples_indices = set()
        for pc_index, std_num in pcstds:
            logging.info("finding samples with std higher than %d stds or lower than %d stds on pc number %d..." % (std_num, -1 * std_num, pc_index))
            pc = pca_out.P[:,pc_index-1] # user start counting from 1 and python from 0
            std_pc = std(pc)
            maxpcstds_samples_indices.update(where((pc > std_num * std_pc) | (pc < -1 * std_num * std_pc))[0])

        if maxpcstds_samples_indices:
            logging.info("excluding samples with max std...")
            self.remove_samples_indices(list(maxpcstds_samples_indices))

    
    def remove_missing_values_sites(self, missing_values_th = 0.03):
        pass
        # was not tested!! 
        # nan are not supported for version 1.0
        # """
        # remove sites that have many missing values
        # many is self.missing_values_th from the values
        # missing_values_th is float between 0 and 1
        # """
        # max_missing_values = int(missing_values_th * self.samples_size)
        # nan_quantity_per_site = isnan(self.data).sum(axis=1) 
        # many_nan_indices = where(nan_quantity_per_site > max_missing_values) # TODO add [0]?
        # logging.debug("Removing %s out of %s sites with more than %s missing values" % (len(many_nan_indices), self.sites_size, max_missing_values))
        # self.exclude_sites_indices(many_nan_indices)
        # logging.debug("%s sites were excluded" % len(many_nan_indices)

    def replace_missing_values_by_mean(self):
        pass
        # was not tested!! 
        # nan are not supported for version 1.0
        # """
        # replaces nan values with the mean of the site
        # """
        # logging.debug("Replacing missing values by site's mean")

        # # mean_per_site = self.get_mean_per_site()
        # masked_data = masked_array(self.data, isnan(self.data)) 
        # mean_per_site = average(masked_data, axis=1)  
        # # TODO is masked_data.mask equal to nan_indices? if so, we don't need to run this "where" line and just use masked_data.mask instead of nan_indices
        # nan_indices = where(masked_data.mask) #TODO add [0] ?                 # find nan values indices
        # self.data[nan_indices] = mean_per_site[nan_indices[0]]    # replace nan values by the mean of each site

    def copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)

    def run():
        pass

    def add_covar_datas(self, covardata_list):
        """
        covardata is the covariates data (without the colums of the sample_ids)
        assumes covardata is in the right format and is sorted by sample_ids as datafile is sorted
        """
        if self.covar is not None:
            covardata_list.insert(0, self.covar)
        self.covar = column_stack(tuple(covardata_list))


class MethylationDataLoader(MethylationData):
    """
    TODO add class doc here
    """
    def __init__(self, datafile, phenofile = None, covarfiles = []):
        data, samples_ids, cpgnames = self._load_and_validate_datafile(datafile)
        sites_size, samples_size = data.shape
        phenotype = self._load_and_validate_phenotype(phenofile, samples_size, samples_ids)
        covar = self._load_and_validate_covar(covarfiles, samples_size, samples_ids)
        super(MethylationDataLoader, self).__init__(data, array(samples_ids), array(cpgnames), phenotype, covar)

    def _load_and_validate_file_of_dimentions(self, datafile, dim):
        """
        validates that the file contains a matrix from dimentions dim
        datafile is not None
        """
        if not isinstance(datafile, file):
            datafile = open(datafile, 'r')

                    
        logging.info("loading file %s..." % datafile.name)

        data = common.load_data_file(datafile.name, dim)
        if data is None:
            common.terminate("the file '%s' is not a %sd matrix" % (datafile.name, dim))

        return data

    def _load_and_validate_datafile(self, datafile):
        """
        returns data (type=float without samples ids and cpgnames) and sample_ids list and cpgnames list
        """
        data = self._load_and_validate_file_of_dimentions(datafile, 2)
        samples_ids = data[0,:][1:]  # extract samples ID
        cpgnames = data[:,0][1:]     # extract methylation sites names
        # remove sample ID and sites names from matrix
        # that kind of assignment will create a copy of O[1:,1:]
        # Note that assignment like self.O = O will not create a copy
        try:
            data = data[1:,1:].astype(float) 
        except ValueError:
            common.terminate("file contains values which are not float")
        # must be called after convertion to float
        logging.info("checking for missing values in datafile file...")
        validate_no_missing_values(data) #TODO remove when missing values are supported

        return data, samples_ids, cpgnames


    def _validate_samples_ids(self, data, samples_size, samples_ids):
        """
        reads matrix from matrix_file_path
        validates that the matrix has number of rows as the number of sample ids
        checks that the sample ids in matrix (the first column) are the same ids as in sample_ids list
        and in the same order
        """
        if samples_size != len(samples_ids):
            common.terminate("the file doesn't include all sample ids %s %s"% (len(data) ,len(samples_ids)))

        matrix_sample_ids = data[:,0]
        if not (samples_ids == matrix_sample_ids).all():
            if len(set(samples_ids)^set(matrix_sample_ids)) != 0:
                common.terminate("sample ids are not identical to the sample ids in data file")
            common.terminate("sample ids are not in the same order as in the datafile") # TODO Elior, should we terminate because sample_ids in files are not in the same order as in datafile?


    def _load_and_validate_samples_info(self, samples_info, samples_size, samples_ids):
        """
        samples_info - path to file containing information about samples (matrix where first column is sample_id)
        samples_info assumed to hold path (not None)
        """
        data = self._load_and_validate_file_of_dimentions(samples_info, 2)
        self._validate_samples_ids(data, samples_size, samples_ids)
        # remove sample IDs from matrix
        # that kind of assignment will create a copy of data[:,1]
        # Note that assignment like self.O = O will not create a copy
        try:
            data = data[:,1:].astype(float) # use only the first phenotype
        except ValueError:
            common.terminate("file contains values which are not float" )
        validate_no_missing_values(data)
        return data

    def _load_and_validate_phenotype(self, phenofile, samples_size, samples_ids):
        """
        returns phenotype data (type=float) without samples ids
        """
        if not phenofile:
            return None
        logging.info("validating phenotype file...")
        pheno = self._load_and_validate_samples_info(phenofile, samples_size, samples_ids)
        if len(pheno[0]) != 1:
            logging.warning("more than one phenotype is not supported. will use only the first phenotype (first column)") # TODO remove when supported
            pheno = pheno[:,1]

        return pheno

    def _load_and_validate_covar(self, covarfiles_list, samples_size, samples_ids):
        """
        concatenate the covariates into one matrix (type=float).
        Make sure all have n rows and the ids of the samples are sorted the same order as the data file
        """
        if not covarfiles_list:
            logging.warning("didn't supply covariates file")
            return None
        logging.info("validating covariates files...")
        
        all_covar = column_stack(tuple([self._load_and_validate_samples_info(covariates, samples_size, samples_ids) for covariates in covarfiles_list]))
        return all_covar

    def add_covar_files(self, covarfiles_list):
        covardata_list = [self._load_and_validate_samples_info(covariates, self.samples_size, self.samples_ids) for covariates in covarfiles_list]
        self.add_covar_datas(covardata_list)

    def upload_new_covaritates_files(self, covarfiles_list):
        """
        overloads self.covar with the covarfiles_list
        """
        self.covar = self._load_and_validate_covar(covarfiles_list, self.samples_size, self.samples_ids)

    def upload_new_phenotype_file(self, phenofile):
        self.phenotype = self._load_and_validate_phenotype(phenofile, self.samples_size, self.samples_ids)

