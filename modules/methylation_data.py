import os
import sys
import copy
import logging
from pickle import dump
from numpy import delete, isnan, where, column_stack, std, array, savetxt, in1d, mean, hstack, vstack
from module import Module
from utils import common, pca, LinearRegression, sitesinfo
from bisect import bisect_right

DEFAULT_PREFIX = "methylation_data"
GLINT_FILE_SUFFIX = "glint" #TODO move to a config file
DATA_SUFFIX = "txt"
DEAFULT_COVAR_NAME = "c"
DEFAULT_PHENO_NAME = "p"
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
    def __init__(self, meth_data, samples_list, sites_list, phenotype = None ,covar = None, covarnames = None, phenonames = None):
        self.data = meth_data
        self.samples_ids = samples_list
        self.cpgnames = sites_list  # the last index of covariate loaded without name, this index wil be added to the default name
        self.sites_size, self.samples_size = self.data.shape
        if (len(self.samples_ids) != self.samples_size):
            common.terminate("got data with %s samples but %s samples ids" % (self.samples_size, len(self.samples_ids)))


        if (len(self.cpgnames) != self.sites_size):
            common.terminate("got data with %s sites but %s cpgnames" % (self.sites_size, len(self.cpgnames)))

        logging.debug("got methylation data with %s sites and %s samples id" % (self.sites_size, self.samples_size))

        self.phenotype = phenotype
        self.covar = covar
        self.covarnames = covarnames
        self.phenonames = phenonames

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
                common.terminate("all sites are about to be remove") 

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
                common.terminate("all samples are about to be remove")

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
        # # the commentes suppose to handle missing values but we dont support NAs now . they calcs the right mean (tested) but it wasnt tested with missing values
        # masked_data = masked_array(self.data, isnan(self.data)) # create masked array
        # return average(masked_data, axis=1)
        return mean(self.data, axis=1)


    def include(self, include_list):
        """
        include_list - a list with sites names to keep in data (remove all other sites)
        this function removes the cpg sites not found in include_list list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("including sites...")
        remove_indices_list = where(False == in1d(self.cpgnames , include_list))[0]
        logging.info("include sites: %s CpGs from the reference list of %s CpGs will be included" % (len(self.cpgnames) - len(remove_indices_list), len(include_list)))
        self.exclude_sites_indices(remove_indices_list)
        logging.debug("methylation data new size is %s sites by %s samples" % self.data.shape)

    def exclude(self, exclude_list):
        """
        exclude_list - a list with sites name to remove from data
        this function removes the cpg sites found in self.exclude list from the data
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        logging.info("excluding sites...")
        indices_list = where(in1d(self.cpgnames , exclude_list))[0]
        self.exclude_sites_indices(indices_list)
        logging.debug("methylation data new size is %s sites by %s samples" % self.data.shape)

    def keep(self, keep_list):
        """
        keep_list- a list of samples ids to keep in data (remove all other sample ids)
        this function removes the samples ids not found in keep_list list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("keeping only samples in the file...")
        remove_indices_list = where(False == in1d(self.samples_ids , keep_list))[0]
        self.remove_samples_indices(remove_indices_list)
        logging.debug("methylation data new size is %s sites by %s samples" % self.data.shape)

    def remove(self, remove_list):
        """
        remove_list - a list of samples ids to remove from data
        this function removes the samples ids found in remove_list from the data
        it updates the samples_size and the samples_ids list 
        """
        logging.info("removing the samples from the file...")
        indices_list = where(in1d(self.samples_ids , remove_list))[0]
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

    def save_sites_and_samples(self, prefix):
        """
        save samples ids withe the covariates and phenotype
        save cpgnames with  information on each site
        """
        sites_filename = prefix + ".sites." + DATA_SUFFIX
        logging.info("Saving cpg names and info to %s" % sites_filename)
        sites_info = sitesinfo.SitesInfoGenerator(self.cpgnames)
        sites_data = column_stack((self.cpgnames, sites_info.chromosomes, sites_info.positions, sites_info.genes, sites_info.categories))
        savetxt(sites_filename, sites_data, delimiter='\t', fmt = '%-12s\t%-4s\t%-12s\t%-22s\t%-22s', header = "cpgid, chromosome, position, gene, category")
        
        samples_filename = prefix + ".samples." + DATA_SUFFIX
        logging.info("Saving samples ids and their phenotype and covariates to %s" % samples_filename)
        samples_data = self.samples_ids
        samples_header = "sampleid"
        fmt = '%-12s'
        
        if self.phenotype is not None:
            samples_data = column_stack((samples_data, self.phenotype))
            pheno_num = self.phenotype.shape[1]
            if pheno_num > 0:
                samples_header += ", " + ", ".join(self.phenonames)
            fmt += '\t%-12s' * pheno_num

        if self.covar is not None:
            samples_data = column_stack((samples_data, self.covar))
            covar_num = self.covar.shape[1]
            if covar_num > 0:
                samples_header += ", " + ", ".join(self.covarnames)
            fmt += '\t%-12s' * covar_num

        savetxt(samples_filename, samples_data, fmt = fmt, header = samples_header, delimiter='\t')

    def save_raw_data(self, prefix = DEFAULT_PREFIX):
        """
        save the data matrix to a text visible file
        save samples ids withe the covariates and phenotype
        save cpgnames with  information on each site
        """
        if prefix is None:
            prefix = DEFAULT_PREFIX
        
        self.save_sites_and_samples(prefix)

        methylation_data_filename = prefix + "." + DATA_SUFFIX


        data = column_stack((self.cpgnames, self.data))
        header = ["ID"]+list(self.samples_ids)

        data = vstack((header, data))

        logging.info("Saving methylation data to %s" % methylation_data_filename)
        savetxt(methylation_data_filename, data, delimiter='\t', fmt='%s')


    def save_serialized_data(self, prefix = DEFAULT_PREFIX):
        """
        serializes this object and saves it to methylation_data_filename
        assumes that methylation_data_filename is a valid file 
        """
        if prefix is None:
            prefix = DEFAULT_PREFIX

        self.save_sites_and_samples(prefix)
        
        filename = prefix + "." + GLINT_FILE_SUFFIX
        with open(filename, 'wb') as f:
            logging.info("Saving methylation data as glint format at %s" % filename)
            dump(self, f)
        f.close()


    def remove_lowest_std_sites(self, lowest_std_th = 0.02):
        """
        input: lowest_std_th threshold for excluding low variance sites, all sites with std lower than lowest_std_th will be excluded 
        lowest_std_th is float between 0 and 1
        """
        logging.info("excluding site with variance lower than %s..." % lowest_std_th)
        # get std for each site
        # sites_std = nanstd(self.data, axis=1) # calc variance consider NaN - missing values are not supported right noe
        sites_std = std(self.data, axis=1) # calc variance - missing values are not allowed

        # sort std and get sites index for each std (sorted, so indices of the lowest std sites will be to the left) 
        std_sorted_indices = sites_std.argsort() # sort the sites_variance and return an array that holds the indices of the sorted values
        # get std list sorted
        std_sorted = sites_std[std_sorted_indices]
        # get the first index in the sorted list which have std higher than lowest_std_th and include all indices started from it
        include_from_index = bisect_right(std_sorted, lowest_std_th)
        if (include_from_index == self.sites_size):
            common.terminate("the provided std parameter excludes all sites (%s)" % lowest_std_th)
        if (include_from_index == 0):
            logging.warning("the provided std parameter excludes no sites (%s)" % lowest_std_th)
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

    
    # def remove_missing_values_sites(self, missing_values_th = 0.03):
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

    # def replace_missing_values_by_mean(self):
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

    def get_phenotypes_indicis(self, names_list):
        if self.phenonames is None: # there are no phenotypes
            common.terminate("there is no phenotype in the data. add phenotype files with --phenofile")
        indices = where(in1d(self.phenonames, names_list))[0]
        if len(indices) != len(names_list):
            common.terminate("some phenotypes names does not match the names in the data")
        return indices

    def get_covariates_indicis(self, covariates_names_list):
        if self.covarnames is None: # there are no covariates
            common.terminate("there is no covariate in the data. add covariate files with --covarfile")
        indices = where(in1d(self.covarnames , covariates_names_list))[0]
        if len(indices) != len(covariates_names_list):
            common.terminate("some covariates names does not match the names in the data")
        return indices

    def get_phenotype_subset(self, names_list):
        if names_list == None:
            logging.info("ignoring phenotype")
            return None
        elif names_list == []:
            logging.info("using all phenotypes")
            return self.phenotype
        else:
            logging.info("using phenotypes %s" % ", ".join(names_list))
            return self.phenotype[:, self.get_phenotypes_indicis(names_list)]
    
    def get_covariates_subset(self, names_list):
        if names_list == None:
            return None
        elif names_list == []:
            logging.info("using all covariates")
            return self.covar
        else:
            logging.info("using covariates %s" % ", ".join(names_list))
            return self.covar[:, self.get_covariates_indicis(names_list)]

    def regress_out(self, matrix_to_regress):
        """
        ** THIS FUNCTION CHANGES THE METHYLATION DATA **

        matrix_to_regress - is n X p where n is number of samples (not None)

        regress out matrix_to_regress and update the data
        """
        residuals =  LinearRegression.regress_out(self.data.transpose(), matrix_to_regress)
        residuals = residuals.transpose()
        self.data = residuals
        
class TitleManager(object):
    def __init__(self):
        self.title_indexes = dict()

    def get_next_index(self, title_base_name):
        if title_base_name not in self.title_indexes:
            self.title_indexes[title_base_name] = 1
        return self.title_indexes[title_base_name]

    def increase_next_index(self, title_base_name, count):
        """
        assumes title_base_name exists
        """
        self.title_indexes[title_base_name] = self.title_indexes[title_base_name] + count

    def generate_title(self, title_base_name, count):
        start_i = self.get_next_index(title_base_name)
        titles_list = ["%s%d" % (title_base_name, i) for i in range(start_i, start_i + count)]
        self.increase_next_index(title_base_name, count)
        return titles_list


class MethylationDataLoader(MethylationData):
    """
    responsible for loading the methylation data files, which means it validates that:
    - the data file:
        - does not contain missing values
        - it is a 2D matrix
        - does not contain values which are not int or float (except the header)
        - Note that it assumes 
             - that there is a header (first row) with the samples names
             - that the first column is cpg names
    - the phenotype file (if supplied):
        - contains the same samples and as the data file
        - the first column is the samples ids (names)
        - Note that only one phenotype is supported - therefore if there is more htan one phenotype only the first will be used
    - the covariates file (is supplied):
        - contains the same samples as the data file

    if one of the condition not fulfilled the program terminates.

    after validation it creates the MethylationData object
    """
    def __init__(self, datafile, phenofile = [], covarfiles = []):
        self.header_manager = TitleManager()
        self.covarnames = None
        self.phenonames = None

        data, samples_ids, cpgnames = self._load_and_validate_datafile(datafile)
        sites_size, samples_size = data.shape
        phenotype, phenonames = self._load_and_validate_phenotype(phenofile, samples_size, samples_ids)
        covar, covarnames = self._load_and_validate_covar(covarfiles, samples_size, samples_ids)
        super(MethylationDataLoader, self).__init__(data, array(samples_ids), array(cpgnames), phenotype, covar, covarnames, phenonames)

    def _load_and_validate_file_of_dimentions(self, datafile, dim):
        """
        validates that the file contains a matrix from dimentions dim
        datafile is not None
        """
        if not isinstance(datafile, file):
            datafile = open(datafile, 'r')
  
        logging.info("loading file %s..." % datafile.name)

        data, samples_ids, cpgnames = common.load_data_file(datafile.name, dim)
        if data is None:
            common.terminate("there is a problem with the format of the file '%s'" % datafile.name)
        if cpgnames is None:
            common.terminate("there are no cpgnames for the sites in the datafile '%s'" % datafile.name)
        if samples_ids is None:
            common.terminate("there are no samples ids (header) in the datafile '%s'" % datafile.name)
        return data, samples_ids, cpgnames

    def _load_and_validate_datafile(self, datafile):
        """
        returns data (type=float without samples ids and cpgnames) and sample_ids list and cpgnames list
        """
        data, samples_ids, cpgnames = self._load_and_validate_file_of_dimentions(datafile, 2)
        logging.info("checking for missing values in datafile file...")
        validate_no_missing_values(data) #TODO remove when missing values are supported

        return data, samples_ids, cpgnames


    def _validate_samples_ids(self, matrix_sample_ids, samples_ids):
        """
        reads matrix from matrix_file_path
        validates that the matrix has number of rows as the number of sample ids
        checks that the sample ids in matrix (the first column) are the same ids as in sample_ids list
        and in the same order
        """
        if not ((samples_ids.size == matrix_sample_ids.size) and ((samples_ids == matrix_sample_ids).all())):
            if len(set(samples_ids)^set(matrix_sample_ids)) != 0:
                common.terminate("sample ids are not identical to the sample ids in data file")
            common.terminate("sample ids are not in the same order as in the datafile") 

    def _load_and_validate_samples_info(self, data_with_samples_info, samples_size, samples_ids):
        """
        data_with_samples_info - path to file containing information about samples (matrix where first column is sample_id)
        data_with_samples_info assumed to hold path (not None)
        """
        data, header, new_samples_ids = self._load_and_validate_file_of_dimentions(data_with_samples_info, 2)
        self._validate_samples_ids(new_samples_ids, samples_ids)
        validate_no_missing_values(data)

        return data, header
    
    def _load_and_validate_samples_data(self, data_list, samples_size, samples_ids, defaule_header_name):
        all_data = []
        all_headers = []

        for datafile in data_list:
            data, header = self._load_and_validate_samples_info(datafile, samples_size, samples_ids)
            if header is None:
                header = self.header_manager.generate_title(defaule_header_name, data.shape[1])
            all_data.append(data)
            all_headers.append(header)

        data = column_stack(tuple(all_data))
        header = hstack(tuple(all_headers))
        return data, header

    def _load_and_validate_phenotype(self, phenofile_list, samples_size, samples_ids, default_pheno_name = DEFAULT_PHENO_NAME):
        """
        returns phenotype data (type=float) without samples ids
        """
        if not phenofile_list:
            return None, None

        logging.info("validating phenotype file...")
        pheno, header = self._load_and_validate_samples_data(phenofile_list, samples_size, samples_ids, default_pheno_name)
        logging.info("new phenotypes found: %s" % ", ".join(header))
        return pheno, header
 

    def _load_and_validate_covar(self, covarfiles_list, samples_size, samples_ids, default_covar_name = DEAFULT_COVAR_NAME):
        """
        concatenate the covariates into one matrix (type=float).
        Make sure all have n rows and the ids of the samples are sorted the same order as the data file
        """
        if not covarfiles_list:
            return None, None

        logging.info("validating phenotype file...")
        covar, header = self._load_and_validate_samples_data(covarfiles_list, samples_size, samples_ids, default_covar_name)
        logging.info("new covariates found: %s" % ", ".join(header))
        
        return covar, header

    def update_covar_header(self, titles_to_add):
        if self.covarnames is not None:
            mutual_names = set(self.covarnames).intersection(set(titles_to_add))
            if mutual_names:
                common.terminate("more than one covariate with the name %s" % str (mutual_names))
            self.covarnames = hstack((self.covarnames, titles_to_add))
        else:
            self.covarnames = titles_to_add
    
    def update_covar_data(self, new_data):
        if self.covar is not None:
            self.covar = column_stack((self.covar, new_data))
        else:
            self.covar = new_data

    def update_pheno_header(self, titles_to_add):
        if self.phenonames is not None:
            mutual_names = set(self.phenonames).intersection(set(titles_to_add))
            if mutual_names:
                common.terminate("more than one phenotype with the name %s" % str (mutual_names))
            self.phenonames = hstack((self.phenonames, titles_to_add))
        else:
            self.phenonames = titles_to_add

    def update_pheno_data(self, new_data):
        if self.phenotype is not None:
            self.phenotype = column_stack((self.phenotype, new_data))
        else:
            self.phenotype = new_data

    def add_covar_datas(self, covardata, default_covar_name = DEAFULT_COVAR_NAME, covarsnames = None):
        """
        adds new covariates to the existing covariates 
        covardata - matrix with the covariates data (without the colums of the sample_ids and the title of the covarnames)
        covarsnames - the names of the covariates in the data. if not supplied a default name will be set
        default_covar_name - a default name for the covariates (in case covarsnames is None)

        assumes covardata is in the right format and the sample_ids order is as in the datafile 
        """
        if covarsnames is None:
            covarsnames = self.header_manager.generate_title(default_covar_name, covardata.shape[1])
        
        self.update_covar_data(covardata)
        self.update_covar_header(covarsnames)
        logging.info("added covariates %s" % ", ".join(covarsnames))

    def add_covar_files(self, covarfiles_list, default_covar_name = DEAFULT_COVAR_NAME):
        """
        adds covariates from a list of covariates files to the glint meth object
        assumes self.covars is not None
        """
        covars, covarsnames = self._load_and_validate_covar(covarfiles_list, self.samples_size, self.samples_ids, default_covar_name)
        self.update_covar_data(covars)
        self.update_covar_header(covarsnames)

    def add_pheno_files(self, phenofiles_list, default_pheno_name = DEFAULT_PHENO_NAME):
        """
        adds phenotypes from a list of phenotype files to the glint meth object
        assumes self.phenotype is not None
        """
        phenos, names = self._load_and_validate_phenotype(phenofiles_list, self.samples_size, self.samples_ids, default_pheno_name)
        self.update_pheno_data(phenos)
        self.update_pheno_header(names)

    def upload_new_covaritates_files(self, covarfiles_list, default_covar_name = DEAFULT_COVAR_NAME):
        """
        overloads self.covar with the covarfiles_list
        """
        self.covar, self.covarnames = self._load_and_validate_covar(covarfiles_list, self.samples_size, self.samples_ids, default_covar_name)

    def upload_new_phenotype_file(self, phenofile):
        self.phenotype, self.phenonames = self._load_and_validate_phenotype(phenofile, self.samples_size, self.samples_ids)


