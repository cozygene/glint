import os
import sys
import copy
import logging
from numpy import loadtxt, delete, isnan, nanvar, where
from numpy.ma import average, masked_array
from module import Module

COMPRESSED_FILENAME = "methylation_data.glint"

class MethylationData( Module ):
    """
    includefile and excludefile cannot be set together
    keepfile and removefile cannot be set together
    TODO add class doc here
    datafile,
                 missing_values_th -  if set, the percentage of sites with the highest number of missing values to remove. a number between 0 and 1
                 lowest_variance_th - if set, the percentage of sites with the lowest variance to remove. a number between 0 and 1
                 min_value -          if set, will remove sites with mean value lower than this. a number between 0 and 1
                 max_value -          if set, will remove sites with mean value higher than this. a number between 0 and 1
                 includefile -        file with a list of sites name to include in the data
                 excludefile -        file with a list of sites name to exclude from data
                 keepfile -           file with list of samples_ids to keep in data
                 removefile -         file with list of samples ids to remove from data
                 glint_data_filename -               the file to save the serialized methylation data, "Glint" format
    """
    def __init__(self, datafile,
                 missing_values_th = 0.3,
                 lowest_variance_th = 0.5,
                 min_value = None,
                 max_value = None,
                 includefile = None,
                 excludefile = None,
                 keepfile = None,
                 removefile = None,
                 glint_data_filename = None):
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

        self.include = self._validate_cpgnames_in_file(includefile)              if includefile else None
        self.exclude = self._validate_cpgnames_in_file(excludefile)              if excludefile else None
        self.keep    = self._validate_sample_ids_in_file(keepfile)               if keepfile    else None
        self.remove  = self._validate_sample_ids_in_file(removefile)             if removefile  else None
        self.missing_values_th  = self._validate_percentage(missing_values_th)   if missing_values_th else None
        self.lowest_variance_th = self._validate_percentage(lowest_variance_th)  if lowest_variance_th else None
        self.glint_data_filename = self._validate_file_path(glint_data_filename) if glint_data_filename else None
        self.min_value = self._validate_methylation_value(min_value)             if min_value else None
        self.max_value = self._validate_methylation_value(max_value)             if max_value else None
        self._validate_min_and_max_mean_values()
       

        # arguments that will hold results of heavy calculations
        # will be good in case a function is called twice
        self._mean_per_site = None

        # run everything
        self.run()

    def _validate_percentage(self, num):
        if not (num <= 1 and num >=0):
            logging.error("must be a percentage (float number between 0 and 1)")
            sys.exit(2) 
        return num

    def _validate_methylation_value(self, num):
        if not (num <= 1 and num >=0):
            logging.error("must be a standard methylation value (float number between 0 and 1)")
            sys.exit(2) 
        return num

    def _validate_min_and_max_mean_values(self):
        """
        validates that the min_value is not greater than max_value
        must be called after self.min_value and self.max_value are set
        """
        if self.min_value is not None and self.max_value is not None:
            if self.max_value <= self.min_value:
                logging.error("min value %s is greater than max value %s" % (self.min_value, self.max_value))
                sys.exit(2) 

   
    def _validate_cpgnames_in_file(self, filepath):
        """
        validates that all cpgs are uniqe and appears in datafile
        """
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

    
    def _validate_sample_ids_in_file(self, filepath):
        """
        validates that all sample ids are uniqe and appears in datafile
        """
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

    def _validate_file_path(self, filepath):
        if not os.path.exists(filepath) :
            logging.error("The file '%s' doesn't exist. Exiting" % filepath)
            sys.exit(2)

    def _validate_file_of_dimentions(self, filepath, dim):
        """
        validates that a file exists and that it is a matrix from dimentions dim
        """
        if filepath is None:
            return None

        self._validate_file_path(filepath)

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
        """
        this function removes the cpg sites not found in self.include list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        if self.include is not None:
            indices_list = [i for i, site in enumerate(self.cpgnames) if site in self.include]
            self.data = self.data[indices_list, :]
            self.cpgnames = self.cpgnames[indices_list]
            self._mean_per_site = self._mean_per_site[indices_list] if self._mean_per_site is not None else None
            self.sites_size = len(self.cpgnames)

    def _exclude_sites_from_data(self, sites_indicies_list):
        """
        this function removes from the data the cpg sites which indices found in sites_indicies_list
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        self.data = delete(self.data, sites_indicies_list, axis = 0)
        self.cpgnames = delete(self.cpgnames, sites_indicies_list)
        self._mean_per_site = delete(self._mean_per_site, indices_list) if self._mean_per_site is not None else None
        self.sites_size = len(self.cpgnames)

    def _exclude(self):
        """
        this function removes the cpg sites found in self.exclude list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        if self.exclude is not None:
            indices_list = [i for i, site in enumerate(self.cpgnames) if site in self.exclude]
            self._exclude_sites_from_data(indices_list)

    def _keep(self):
        """
        this function removes the samples ids not found in self.keep list from the data
        it updates the samples_size and the samples_ids list 
        """
        if self.keep is not None:
            indices_list = [i for i, id in enumerate(self.samples_ids) if site in self.keep]
            self.data =  self.data[:, indices_list] 
            self.samples_ids = self.samples_ids[indices_list]
            self.samples_size = len(self.samples_ids)

    def _remove(self):
        """
        this function removes the samples ids found in self.remove list from the data
        it updates the samples_size and the samples_ids list 
        """
        if self.remove is not None:
            indices_list = [i for i, id in enumerate(self.samples_ids) if site in self.remove]
            self.data = delete(self.data, indices_list, axis = 1)
            self.samples_ids = delete(self.samples_ids, indices_list)
            self.samples_size = len(self.samples_ids)
    
    def _exclude_sites_with_low_mean(self):
        """
        removes sites with mean < min_value
        """
        if self.min_value is not None:
            min_values_indices = where(self.get_mean_per_site() < self.min_value)   
            self._exclude_sites_from_data(min_values_indices)

    def _exclude_sites_with_high_mean(self):
        """
        removes sites with mean > max_value
        """
        if self.max_value is not None:
            max_values_indices = where(self.get_mean_per_site() > self.max_value)   
            self._exclude_sites_from_data(max_values_indices)

    def _save(self):
        if self.glint_data_filename is not None:
            import pickle
            with open(self.glint_data_filename, 'wb') as f:
                pickle.dump(self, f)


    def run(self):
        self._include()
        self._exclude()
        self._keep()
        self._remove()
        self._exclude_sites_with_high_mean()
        self._exclude_sites_with_low_mean()
        self._save()

    def remove_lowest_variance_sites(self):
        """
        removes the sites with the lowest variance.
        removes self.lowest_variance_th from the sites.
        """
        if self.lowest_variance_th is not None:
            sites_variance = nanvar(self.data, axis=1) # calc variance consider NaN
            var_sorted_indices = sites_variance.argsort() # sort the sites_variance and return an array that holds the indices of the sorted values
            quantity_to_remove = self.lowest_variance_th * self.sites_size
            lowest_variance_indices = var_sorted_indices[:quantity_to_remove]

            self._exclude_sites_from_data(lowest_variance_indices)


    def get_mean_per_site(self):
        """
        returns array that contains the mean value (average) for each methylation site
        if the function called more than once, it won't calculate the average again
        """
        if self._mean_per_site is None: # if this function is called the first time
            masked_data = masked_array(self.data,isnan(self.data)) # create masked array
            self._mean_per_site = average(masked_data, axis=1)     # calc average for each site and save it 

        return self._mean_per_site

    def remove_missing_values_sites(self):
        """
        remove sites that have many missing values
        many is self.missing_values_th from the values
        """
        if self.missing_values_th is not None:
            max_missing_values = self.missing_values_th * self.samples_size
            nan_quantity_per_site = isnan(self.data).sum(axis=1)
            many_nan_indices = where(nan_quantity_per_site > max_missing_values)
            self._exclude_sites_from_data(many_nan_indices)

    def replace_missing_values_by_mean(self):
        """
        replaces nan values with the mean of the site
        """
        nan_indices = where(masked_data.mask)                  # find nan values indices
        self.data[nan_indices] = self.get_mean_per_site()[nan_indices[0]]    # replace nan values by the mean of each site


    def copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)
