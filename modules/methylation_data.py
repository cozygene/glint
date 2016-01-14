import os
import sys
import copy
import logging
from pickle import dump
from numpy import loadtxt, delete, isnan, nanvar, where
from numpy.ma import average, masked_array
from module import Module
from utils import common

COMPRESSED_FILENAME = "methylation_data"

class MethylationData( Module ):
    """
    TODO add class doc here
    """
    def __init__(self, datafile):
        data = self._load_and_validate_file_of_dimentions(datafile, 2) 
        self.samples_ids = data[0,:][1:]                 # extract samples ID
        self.cpgnames = data[:,0][1:]                    # extract methylation sites names

        # remove sample ID and sites names from matrix
        # that kind of assignment will create a copy of O[1:,1:], TODO: do we need a new copy here? i don't think so
        # Note that assignment like self.O = O will not create a copy
        self.data = data[1:,1:].astype(float) 
        self.sites_size, self.samples_size = self.data.shape


        logging.debug("Got methylation data with %s sites and %s samples id" % (self.sites_size, self.samples_size))


    def _validate_file_path(self, filepath):
        if not os.path.exists(filepath) :
            logging.error("The file '%s' doesn't exist. Exiting" % filepath)
            common.terminate(self.__class__.__name__)

    def _load_and_validate_file_of_dimentions(self, filepath, dim):
        """
        validates that a file exists and that it is a matrix from dimentions dim
        """
        if filepath is None:
            return None

        self._validate_file_path(filepath)
        logging.info("Loading file %s..." % filepath)
        data = loadtxt(filepath, dtype = str)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if len(data.shape) != dim:
            logging.error("The file '%s' is not a %sd matrix" % (filepath, dim))
            common.terminate(self.__class__.__name__)

        return data

    def _exclude_sites_from_data(self, sites_indicies_list):
        """
        this function removes from the data the cpg sites which indices found in sites_indicies_list
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        self.data = delete(self.data, sites_indicies_list, axis = 0)
        self.cpgnames = delete(self.cpgnames, sites_indicies_list)
        self.sites_size = len(self.cpgnames)
        logging.debug("methylation data new size is %s" % str(self.data.shape))
        if (self.data.shape[0] != self.sites_size):
            logging.error("After excluding sites, methylation data sites size is %s but we got %s" % (self.data.shape[0], self.sites_size))

    def get_mean_per_site(self):
        """
        returns array that contains the mean value (average) for each methylation site
        if the function called more than once, it won't calculate the average again
        """
        masked_data = masked_array(self.data, isnan(self.data)) # create masked array
        return average(masked_data, axis=1)

    def include(self, include_list):
        """
        this function removes the cpg sites not found in include_list list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("Including %s out of %s sites..." % (len(include_list), self.sites_size))
        indices_list = [i for i, site in enumerate(self.cpgnames) if site in include_list]
        self.data = self.data[indices_list, :]
        self.cpgnames = self.cpgnames[indices_list]
        self.sites_size = len(self.cpgnames)
        logging.debug("methylation data new size is %s" % self.data.shape)
        if (self.data.shape[0] != self.sites_size):
            logging.error("After including sites, methylation data sites size is %s but we got %s" % (self.data.shape[0], self.sites_size))

    def exclude(self, exclude_list):
        """
        this function removes the cpg sites found in self.exclude list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("Excluding %s out of %s sites..." % (len(exclude_list), self.sites_size))
        indices_list = [i for i, site in enumerate(self.cpgnames) if site in exclude_list]
        self._exclude_sites_from_data(indices_list)

    def keep(self, keep_list):
        """
        this function removes the samples ids not found in keep_list list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("Keeping %s out of %s samples..." % (len(keep_list), self.samples_size))
        indices_list = [i for i, id in enumerate(self.samples_ids) if site in keep_list]
        self.data =  self.data[:, indices_list] 
        self.samples_ids = self.samples_ids[indices_list]
        self.samples_size = len(self.samples_ids)
        if (self.data.shape[1] != self.samples_size):
            logging.error("After kepping samples, methylation data samples size is %s but we got %s" % (self.data.shape[1], self.samples_size))

    def remove(self, remove_list):
        """
        this function removes the samples ids found in remove_list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("Removing %s out of %s samples..." % (len(remove_list), self.samples_size))
        indices_list = [i for i, id in enumerate(self.samples_ids) if site in remove_list]
        self.data = delete(self.data, indices_list, axis = 1)
        self.samples_ids = delete(self.samples_ids, indices_list)
        self.samples_size = len(self.samples_ids)
        if (self.data.shape[1] != self.samples_size):
            logging.error("After removing samples, methylation data samples size is %s but we got %s" % (self.data.shape[1], self.samples_size))

    def exclude_sites_with_low_mean(self, min_value):
        """
        removes sites with mean < min_value
        """
        min_values_indices = where(self.get_mean_per_site() < min_value)  
        logging.info("Removing %s out of %s sites with mean lower than %s" % (len(min_values_indices),  self.sites_size, min_value))
        self._exclude_sites_from_data(min_values_indices)

    def exclude_sites_with_high_mean(self, max_value):
        """
        removes sites with mean > max_value
        """
        max_values_indices = where(self.get_mean_per_site() > max_value)
        logging.info("Removing %s out of %s sites with mean greater than %s" % (len(max_values_indices),  self.sites_size, max_value))
        self._exclude_sites_from_data(max_values_indices)

    def save(self, methylation_data_filename):
        """
        serializes this object and saves it to methylation_data_filename
        assumes that methylation_data_filename is a valid file 
        """
        with open(methylation_data_filename, 'wb') as f:
            logging.info("Saving methylation data as glint format at %s" % methylation_data_filename)
            dump(self, f)


    def remove_lowest_variance_sites(self, lowest_variance_th):
        """
        removes the sites with the lowest variance.
        removes self.lowest_variance_th from the sites.
        lowest_variance_th is float between 0 and 1
        """
        sites_variance = nanvar(self.data, axis=1) # calc variance consider NaN
        var_sorted_indices = sites_variance.argsort() # sort the sites_variance and return an array that holds the indices of the sorted values
        quantity_to_remove = int(lowest_variance_th * self.sites_size)
        logging.debug("Removing %s out of %s sites with the lowest variance" % (quantity_to_remove, self.sites_size))
        lowest_variance_indices = var_sorted_indices[:quantity_to_remove]
        self._exclude_sites_from_data(lowest_variance_indices)
    
    def remove_missing_values_sites(self, missing_values_th):
        """
        remove sites that have many missing values
        many is self.missing_values_th from the values
        missing_values_th is float between 0 and 1
        """
        max_missing_values = int(missing_values_th * self.samples_size)
        nan_quantity_per_site = isnan(self.data).sum(axis=1) 
        many_nan_indices = where(nan_quantity_per_site > max_missing_values)
        logging.debug("Removing %s out of %s sites with more than %s missing values" % (len(many_nan_indices), self.sites_size, max_missing_values))
        self._exclude_sites_from_data(many_nan_indices)

    def replace_missing_values_by_mean(self):
        """
        replaces nan values with the mean of the site
        """
        logging.debug("Replacing missing values by site's mean")
        mean_per_site = self.get_mean_per_site()

        # TODO is masked_data.mask equal to nan_indices? if so, we don't need to run this "where" line and just use masked_data.mask instead of nan_indices
        nan_indices = where(masked_data.mask)                    # find nan values indices

        self.data[nan_indices] = mean_per_site[nan_indices[0]]    # replace nan values by the mean of each site


    def copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)

    def run():
        pass
