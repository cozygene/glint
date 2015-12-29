import os
import sys
import logging
from numpy import loadtxt, delete, isnan, nanvar, where
from numpy.ma import average, masked_array
from module import Module

class MethylationData( Module ):
    """
    includefile and excludefile cannot be set together
    keepfile and removefile cannot be set together
    TODO add class doc here
    """
    def __init__(self, datafile, missing_values_th = 0.3, lowest_variance_th = 0.5, includefile = None, excludefile = None, keepfile = None, removefile = None):
        print includefile
        data = self._validate_str_matrix_file(datafile) 
        self.samples_ids = data[0,:][1:]         # extract samples ID
        self.cpgnames = data[:,0][1:]           # extract methylation sites names

        # remove sample ID and sites names from matrix
        # that kind of assignment will create a copy of O[1:,1:], TODO: do we need a new copy here? i don't think so
        # Note that assignment like self.O = O will not create a copy
        self.data = data[1:,1:].astype(float) 
        self.sites_size, self.samples_size = self.data.shape

        if excludefile is not None and includefile is not None:
            logging.error("includefile and excludefile cannot be set together")
            sys.exit(2)

        if keepfile is not None and removefile is not None:
            logging.error("removefile and keepfile cannot be set together")
            sys.exit(2)

        self.include = self._validate_cpgnames_in_file(includefile)  if includefile else None
        self.exclude = self._validate_cpgnames_in_file(excludefile)  if excludefile else None
        self.keep    = self._validate_sample_ids_in_file(keepfile)   if keepfile    else None
        self.remove  = self._validate_sample_ids_in_file(removefile) if removefile  else None

        self.missing_values_th  = self._validate_percentage(missing_values_th)
        self.lowest_variance_th = self._validate_percentage(lowest_variance_th)

        self.run()

    def _validate_percentage(self, num):
        if not (num <= 1 and num >=0):
            print("must be a float between 0 and 1")
            sys.exit(2) 
        return num

    def get_copy_of_data(self):
        """
        returns a copy of the data.
        if you want a copy of the object just run copy.deepcopy(you_object)
        """
        return self.data.copy()
        
    """
    validates that all cpgs are uniqe and appears in datafile
    """
    def _validate_cpgnames_in_file(self, filepath):
        data = self._validate_str_vector_file(filepath)
        data_set = set(data)
        if len(data) != len(data_set):
            logging.error("list of cpgnames in file %s cannot include a cpgname more than once")
            sys.exit(2)

        diff =  data_set.difference(set(self.cpgnames))
        if diff != set([]):
            logging.error("list of cpgnames in file %s includes a cpgnames that are not specified in datafile: %s" % diff)
            sys.exit(2)

        return data

    """
    validates that all sample ids are uniqe and appears in datafile
    """
    def _validate_sample_ids_in_file(self, filepath):
        data = self._validate_str_vector_file(filepath)
        data_set = set(data)
        if len(data) != len(data_set):
            logging.error("list of sample ids in file %s cannot include a cpgname more than once")
            sys.exit(2)

        diff =  data_set.difference(set(self.samples_ids))
        if diff != set([]):
            logging.error("list of sample ids in file %s includes a sample ids that are not specified in datafile: %s" % diff)
            sys.exit(2)

        return data

    """
    validates that a file exists and that it is a matrix from dimentions dim
    """
    def _validate_file_of_dimentions(self, filepath, dim):
        if filepath is None:
            return None

        if not os.path.exists(filepath) :
            logging.error("The file '%s' doesn't exist. Exiting" % filepath)
            sys.exit(2)

        data = loadtxt(filepath, dtype = str)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if len(data.shape) != dim:
            logging.error("The file '%s' is not a %sd vector" % (filepath, dim))
            sys.exit(2)

        return data

    def _validate_str_matrix_file(self, filepath):
        return self._validate_file_of_dimentions(filepath, 2)

    def _validate_str_vector_file(self, filepath):
        return self._validate_file_of_dimentions(filepath, 1)

    def _include(self):
        if self.include is not None:
            indices_list = [i for i, site in enumerate(self.cpgnames) if site in self.include]
            self.data = self.data[indices_list, :]
            self.cpgnames = self.cpgnames[indices_list]
            self.sites_size = len(self.cpgnames)

    def _exclude_sites_from_data(self, sites_indicies_list):
        self.data = delete(self.data, sites_indicies_list, axis = 0)
        self.cpgnames = delete(self.cpgnames, sites_indicies_list)
        self.sites_size = len(self.cpgnames)

    def _exclude(self):
        if self.exclude is not None:
            indices_list = [i for i, site in enumerate(self.cpgnames) if site in self.exclude]
            self._exclude_sites_from_data(indices_list)

    def _keep(self):
        if self.keep is not None:
            indices_list = [i for i, id in enumerate(self.samples_ids) if site in self.keep]
            self.data =  self.data[:, indices_list] 
            self.samples_ids = self.samples_ids[indices_list]
            self.samples_size = len(self.samples_ids)

    def _remove(self):
        if self.remove is not None:
            indices_list = [i for i, id in enumerate(self.samples_ids) if site in self.remove]
            self.data = delete(self.data, indices_list, axis = 1)
            self.samples_ids = delete(self.samples_ids, indices_list)
            self.samples_size = len(self.samples_ids)
        
    def run(self):
        self._include()
        self._exclude()
        self._keep()
        self._remove()

    def remove_lowest_variance_sites(self):
        """
        removes the sites with the lowest variance.
        removes self.lowest_variance_th from the sites.
        """
        sites_variance = nanvar(self.data, axis=1) # calc variance consider NaN
        var_sorted_indices = sites_variance.argsort() # sort the sites_variance and return an array that holds the indices of the sorted values
        quantity_to_remove = self.lowest_variance_th * self.sites_size
        lowest_variance_indices = var_sorted_indices[:quantity_to_remove]

        self._exclude_sites_from_data(lowest_variance_indices)

    def remove_missing_values_sites(self):
        """
        remove sites that have many missing values
        many is self.missing_values_th from the values
        """

        self.missing_values_th = 0.012 #TODO remove this
        max_missing_values = self.missing_values_th * self.samples_size
        nan_quantity_per_site = isnan(self.data).sum(axis=1)
        many_nan_indices = where(nan_quantity_per_site > max_missing_values)
        self._exclude_sites_from_data(many_nan_indices)

    def replace_missing_values_by_mean(self):
        """
        replaces nan values with the mean of the site
        """
        masked_data = masked_array(self.data,isnan(self.data))
        sites_mean = average(masked_data, axis=1)
        nan_indices = where(masked_data.mask)
        self.data[nan_indices] = sites_mean[nan_indices[0]]
        print isnan(self.data).sum()

