import os
import sys
import copy
import logging
from pickle import dump
from numpy import loadtxt, delete, isnan, nanstd, where
from numpy.ma import average, masked_array
from module import Module
from utils import common, search
from bisect import bisect_right

COMPRESSED_FILENAME = "methylation_data"

def validate_no_missing_values(data):
    """
    nan are not supported for version 1.0
    Note: data must be number and not string
    """
    if  isnan(data).sum() > 0:
        common.terminate("missing values are not supported at this version")


class MethylationData( Module ):
    """
    TODO add class doc here
    """
    def __init__(self, datafile, phenofile = None, covarfile = None):

        self.data, self.samples_ids, self.cpgnames = self._load_and_validate_datafile(datafile)        
        self.sites_size, self.samples_size = self.data.shape
        logging.debug("got methylation data with %s sites and %s samples id" % (self.sites_size, self.samples_size))

        self.phenotype = self._load_and_validate_phenotype(phenofile)
        self.covar = self._load_and_validate_covar(covarfile)

    def _load_and_validate_file_of_dimentions(self, datafile, dim):
        """
        validates that the file contains a matrix from dimentions dim
        datafile is not None
        """
        if not isinstance(datafile, file):
            datafile = open(datafile, 'r')
        logging.info("loading file %s..." % datafile.name)

        data = loadtxt(datafile, dtype = str)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if data.ndim != dim:
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
        data = data[1:,1:].astype(float) 
        
        # must be called after convertion to float
        logging.info("checking for missing values in datafile file...")
        validate_no_missing_values(data) #TODO remove when missing values are supported

        return data, samples_ids, cpgnames

    """
    reads matrix from matrix_file_path
    validates that the matrix has number of rows as the number of sample ids
    checks that the sample ids in matrix (the first column) are the same ids as in sample_ids list
    and in the same order
    """
    def _validate_samples_ids(self, data):
        if len(data) != len(self.samples_ids):
            common.terminate("the file doesn't include all sample ids")

        matrix_sample_ids = data[:,0]

        if not (self.samples_ids == matrix_sample_ids).all(): # TODO check this is not by order
            if len(set(self.samples_ids)^set(matrix_sample_ids)) != 0:
                common.terminate("sample ids are not identical to the sample ids in data file")
            common.terminate("sample ids are not in the same order as in the datafile") # TODO should we terminate because of order?

    def _load_and_validate_samples_info(self, samples_info):
        """
        samples_info - path to file containing information about samples (matrix where first column is sample_id)
        samples_info assumed to hold path (not None)
        """
        data = self._load_and_validate_file_of_dimentions(samples_info, 2)
        self._validate_samples_ids(data)

        # remove sample IDs from matrix
        # that kind of assignment will create a copy of data[:,1]
        # Note that assignment like self.O = O will not create a copy
        data = data[:,1:].astype(float) # use only the first phenotype # TODO should check if can convert  to float
        validate_no_missing_values(data)
        return data

    def _load_and_validate_phenotype(self, phenofile):
        """
        returns phenotype data (type=float) without samples ids
        """
        if not phenofile:
            return None
        logging.info("validating phenotype file...")
        pheno = self._load_and_validate_samples_info(phenofile)
        if len(pheno[0]) != 1:
            logging.warning("more than one phenotype is not supported. will use only the first phenotype (first column)") # TODO remove when supported
            pheno = pheno[:,1] 
        return pheno

    def _load_and_validate_covar(self, covariates):
        """
        returns covariates data (type=float) without samples ids
        """
        if not covariates:
            logging.warning("didn't supply covariates file")
            return None
        logging.info("validating covariates file...")
        return self._load_and_validate_samples_info(covariates)

    def _exclude_sites_from_data(self, sites_indicies_list):
        """
        this function removes from the data the cpg sites which indices found in sites_indicies_list
        it updates the sites_size, the cpgnames list
        """
        self.data = delete(self.data, sites_indicies_list, axis = 0)
        self.cpgnames = delete(self.cpgnames, sites_indicies_list)
        self.sites_size = len(self.cpgnames)
        if (self.data.shape[0] != self.sites_size): # TODO remove this after tests
            common.terminate("After excluding sites, methylation data sites size is %s but we got %s" % (self.data.shape[0], self.sites_size))

    def _exclude_samples_from_data(self, indices_list):
        """
        removes from the data the samples which indices found in indices_list
        it updates the samples_size, the sample_ids lists
        it also remove the sample from phenotype and covariates file if provided
        """
        self.data = delete(self.data, indices_list, axis = 1)
        if self.phenotype is not None:
            self.phenotype = delete(self.phenotype, indices_list, axis = 0)
        if self.covar is not None:
            self.covar = delete(self.covar, indices_list, axis = 0)
        self.samples_ids = delete(self.samples_ids, indices_list)
        self.samples_size = len(self.samples_ids)
        if (self.data.shape[1] != self.samples_size): # TODO remove this after tests
            common.terminate("After removing samples, methylation data samples size is %s but we got %s" % (self.data.shape[1], self.samples_size))
    
    def get_mean_per_site(self):
        """
        returns array that contains the mean value (average) for each methylation site
        """
        masked_data = masked_array(self.data, isnan(self.data)) # create masked array
        return average(masked_data, axis=1)

    def include(self, include_list):
        """
        this function removes the cpg sites not found in include_list list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("including sites...")
        indices_list = [i for i, site in enumerate(self.cpgnames) if site in include_list]
        self.data = self.data[indices_list, :]
        # TODO remove this test if it didnt fail
        if sorted(indices_list) != indices_list:
            common.terminate("indices_list must be sorted in order to do self.cpgnames[indices_list]")
        self.cpgnames = self.cpgnames[indices_list]
        self.sites_size = len(self.cpgnames)
        logging.debug("methylation data new size is %s" % str(self.data.shape))
        if (self.data.shape[0] != self.sites_size):
            common.terminate("After including sites, methylation data sites size is %s but we got %s" % (self.data.shape[0], self.sites_size))
        logging.debug("%s sites were included" % len(indices_list))

    def exclude(self, exclude_list):
        """
        this function removes the cpg sites found in self.exclude list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("excluding sites...")
        indices_list = [i for i, site in enumerate(self.cpgnames) if site in exclude_list]
        self._exclude_sites_from_data(indices_list)
        logging.debug("%s sites were excluded" % len(indices_list))

    def keep(self, keep_list):
        """
        this function removes the samples ids not found in keep_list list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("keeping samples...")
        remove_indices_list = [i for i, site in enumerate(self.samples_ids) if site not in keep_list]
        self._exclude_samples_from_data(remove_indices_list)
        logging.debug("%s samples were removed" % len(remove_indices_list))

    def remove(self, remove_list):
        """
        this function removes the samples ids found in remove_list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("removing samples...")
        indices_list = [i for i, site in enumerate(self.samples_ids) if site in remove_list]
        self._exclude_samples_from_data(indices_list)
        logging.debug("%s samples were removed" % len(remove_list))

    def exclude_sites_with_low_mean(self, min_value):
        """
        removes sites with mean < min_value
        """
        logging.info("excluding sites with mean lower than %s..." % min_value)
        min_values_indices = where(self.get_mean_per_site() < min_value)  
        self._exclude_sites_from_data(min_values_indices)
        logging.debug("%s sites were excluded" % len(min_values_indices))

    def exclude_sites_with_high_mean(self, max_value):
        """
        removes sites with mean > max_value
        """
        logging.info("removing sites with mean greater than %s..." % max_value)
        max_values_indices = where(self.get_mean_per_site() > max_value)
        self._exclude_sites_from_data(max_values_indices)
        logging.debug("%s sites were excluded" % len(max_values_indices))

    def save(self, methylation_data_filename):
        """
        serializes this object and saves it to methylation_data_filename
        assumes that methylation_data_filename is a valid file 
        """
        with open(methylation_data_filename, 'wb') as f:
            logging.info("Saving methylation data as glint format at %s" % methylation_data_filename)
            dump(self, f)

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
        self._exclude_sites_from_data(std_sorted_indices[:include_from_index])
        logging.debug("%s sites were excluded" % len(std_sorted_indices[:include_from_index]))


    
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
        # many_nan_indices = where(nan_quantity_per_site > max_missing_values)
        # logging.debug("Removing %s out of %s sites with more than %s missing values" % (len(many_nan_indices), self.sites_size, max_missing_values))
        # self._exclude_sites_from_data(many_nan_indices)
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
        # nan_indices = where(masked_data.mask)                    # find nan values indices
        # self.data[nan_indices] = mean_per_site[nan_indices[0]]    # replace nan values by the mean of each site

    def copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)

    def run():
        pass
