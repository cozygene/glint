import logging
from pickle import load
from modules import methylation_data
from utils import common
import argparse
from module_parser import ModuleParser
from numpy import loadtxt
import os


# cpgs on X and Y chromosomes, which increase the gender signal in refactor
HUMAN_X_Y = os.path.join(os.path.dirname(__file__),"assets/HumanMethylationSites_X_Y.txt")
# list of nonspecific methylation sites - cpgs that are caught by nonspecific phrobes (can catch few cpgs)
NONSPECIFIC_PROBES = os.path.join(os.path.dirname(__file__),"assets/nonspecific_probes.txt")
# list of methylation sites (cpgs) that are also SNPs. connecting between the genetics (and therfore ancestry) and the methylation
POLYMORPHIC_CPGS = os.path.join(os.path.dirname(__file__),"assets/polymorphic_cpgs.txt")     

class MethylationDataParser(ModuleParser): 

    def __init__(self, parser):
        """
        --datafile: is a path to a data file which can be matrix or .glint file.
        --phenofiles: list of the phenotypes files
        --covarfiles: list of covariates files. It can get more than one covariates file by supplying a list of files:
                        glint.py --covar path_to_covar_file1 path_to_covar_file2 ...
                    Note that if you supply the flag more than once:
                        glint.py --covar path_to_covar_file1 --covar path_to_covar_file2 ...
                    only file2 will be used as covariates
                    The covariates file should have header specifying each covariate name or ID. If no header is supplied - default value is given as names.
                    If glint file is supplied as datafile (and not text data file), the covariates will be added to it just for the scipt execution. if --gsave is also supplied than the new covariates will be added and saved to the old covariates
        --covarfiles:    list of the covariates names to use. the list is subset of the names of the covariates (the names appear in the header of the covariates file or the default names)
        --maxpcstd: exclude samples with low of high std:
                    receives list of tuples where the first index in each tuple is the pc index and the second index is the std number:
                        glint.py --maxpcstd (index_1,std_1) (index_2,std_2), ...(index_n,std_n)
                    then:
                    for each i from 1 to n:
                        calculates the std of pc index_i, call it S
                        exclude samples with std higher than  std_i*S or lower than -1*std_i*S
        --include:  a path to file with a list of sites to include in the data; removes the rest of the sites.
                    the program validates that there is no site in the list which doesn't appear in the datafile.
        --exclude:  a path to file with a list of sites to exclude from the data; includes the rest of the sites.
                    the program validates that there is no site in the list which doesn't appear in the datafile.
        --keep:     a path to file with a list of samples to include in the data; removes the rest of the samples.
                    the program validates that there is no sample in the list which doesn't appear in the datafile.
        --remove:   a path to file with a list of samples to exclude in the data; includes the rest of the samples.
                    the program validates that there is no sample in the list which doesn't appear in the datafile.
        --minmean:  will remove any site with mean < minmean (terminates if minmean is not a a standard methylation value  between 0 and 1)
        --maxmean:  will remove any site with mean > maxmean (terminates if maxmean is not a a standard methylation value: between 0 and 1)
        --gsave:    save the data.
                    output 3 file:
                        1. a serialized file so next time it will be fater to load it. 
                        2. a list of the samples that are found in the data, also saves the covariates and phenotype of the samples 
                        3. a list of the sites that are found in the data, also saves information about each site: it's chromosome, positions genes and categories.
                    Note that --out flag can be supplied to add a different prefix for each output file
        
        Notes:
            - missing values are not supported for now
            - the program validates that the samples in covariates and phenotype file match the samples in the datafile
            - if a glint file is loaded (with --datafile) and a new covariates or phenotype files are provided (with --covarfiles / --phenofiles) - the new covariates
               and phenotype will be added to those found in the glint. the glint file remains the same unless --gsave is selected which will
               save a new glint file (Note that --out flag can be supplied to add a different prefix for each output file).
            - there are options avaliable for the programmer:
                - get a copy of the methylation data object so if you change the data in one copy the data of the other one remains the same. (function copy)
                - remove samples or sites by list of their indices (functions exclude_sites_indices remove_samples_indices)
            - the program wil terminate if:
                - there are missing values in the data file
                - the provided stdth parameter excludes all sites ( will warn if it excluded no site)
                - if all the sites / samples will be removed after --exclude  / --remove
                - the covariates/phenotype file does not contain all samples that are in the datafile, or contains more samples than in the datafile
        """
        required = parser.add_argument_group('1.Required arguments') # numbering in the group name because help print it by abc order
        required.add_argument('--datafile', type = argparse.FileType('r'), required = True,  help = "A data matrix file of beta-normalized methylation levels or a .glint file")

        optional = parser.add_argument_group('2.Data management options ')
        optional.add_argument('--covarfile', type = argparse.FileType('r'), nargs='*', help = "A covariates file (or files)")
        optional.add_argument('--phenofile', type = argparse.FileType('r'), nargs='*', help = "A phenotype file (or files)")
        optional.add_argument('--maxpcstd', metavar=('PC_INDEX (TODO Elior, change those names?', 'STD_COUNT'), type = int, action = 'append', nargs = 2, help = "TODO Elior, edit: pc index and std number of times for removing outliers")
            
        group1 = optional.add_mutually_exclusive_group(required = False)
        group1.add_argument('--include', type = argparse.FileType('r'), help = "A file with a list of sites to include in the data; removes the rest of the sites")
        group1.add_argument('--exclude', type = argparse.FileType('r'), help = "A file with a list of sites to exclude from the data; includes the rest of the sites")

        group2 = optional.add_mutually_exclusive_group(required = False)
        group2.add_argument('--keep',   type = argparse.FileType('r'), help = "A file with a list of samples to include in the data; removes the rest of the samples")
        group2.add_argument('--remove', type = argparse.FileType('r'), help = "A file with a list of samples to exclude in the data; includes the rest of the samples")

        optional.add_argument('--gsave', action='store_true', help = "Save the data in a glint format; makes following executions faster")
        optional.add_argument('--txtsave', action='store_true', help = "Save the data in a visible format (text)")

        def methylation_value(num):
            num = float(num)
            if not (num <= 1 and num >=0):
                common.terminate("minmean/maxmean must be a standard methylation value (float number between 0 and 1)")
            return num

        optional.add_argument('--minmean', type = methylation_value, help = "A threshold for the minimal mean methylation level to consider")
        optional.add_argument('--maxmean', type = methylation_value, help = "A threshold for the maximal mean methylation level to consider")
        
        optional.add_argument('--rmxy', action='store_true', help = "remove methylation sites on X and Y chromosomes")
        optional.add_argument('--rmns', action='store_true', help = "remove nonspecific sites (nonspecific probes)")
        optional.add_argument('--rmpoly', action='store_true', help = "remove polymorphic sites (CpGs which are also SNPs)")
        
        def std_value(num):
            try:
                num = float(num)
            except:
                common.terminate("minstd must be a float between 0 and 1")
            
            if not (num <= 1 and num >=0):
                common.terminate("minstd must be a float between 0 and 1")
            return num
        optional.add_argument('--minstd',  type = std_value, help = "threshold for excluding low variance sites (all sites with std lower than this threshold will be excluded)")
      
        super(MethylationDataParser, self).__init__(required, optional)
        
    def validate_args(self, args):
        super(MethylationDataParser, self).validate_args(args)
        self._validate_min_and_max_mean_values(args.minmean, args.maxmean)

    def _load_and_validate_ids_in_file(self, fileobj, optional_ids_list, dim = 1):
        """
        validates that the file contains a vector of dimentions dim
        loads a vector file contianing ids list
        warns if there are duplicate ids in the file or if there are ids which are not found in optional_ids_list
        fails if this file is not a list
        """
        if not isinstance(fileobj, file):
            fileobj = open(fileobj, 'r')
        
        logging.info("loading file %s..." % fileobj.name)
        try:
            data = loadtxt(fileobj.name, dtype = str)
        except:
            common.terminate("There was error reading the file '%s', make sure you seperate the values with space, tab or comma" % (fileobj.name, dim))
        
        if data.ndim == 0: # file contains only one item
            data = [data.item()]
        elif data.ndim == 2:
            common.terminate("The file '%s' is not a %sd vector" % (fileobj.name, dim))

        data_set = set(data)
        if len(data) != len(data_set):
            logging.warning("The file %s contains ids more than once" % fileobj.name)

        diff =  data_set.difference(set(optional_ids_list))
        if diff != set([]):
            logging.warning("The file %s contains ids that are not found in the datafile: %s" % (fileobj.name, diff))

        return data

    def _validate_min_and_max_mean_values(self, min_value, max_value):
        """
        validates that the min_value is not greater than max_value
        must be called after self.min_value and self.max_value are set
        """
        if min_value is not None and max_value is not None:
            if max_value <= min_value:
                common.terminate("min value %s is greater than max value %s" % (min_value, max_value))
    
    # must  be called after init_data
    def preprocess_sites_data(self):
        if self.args.include is not None:
            self.module.include(self.include_list)
        if self.args.exclude is not None:
            self.module.exclude(self.exclude_list)
        # exclude min/max values
        if self.args.minmean is not None:
            self.module.exclude_sites_with_low_mean(self.args.minmean)
        if self.args.maxmean is not None:
            self.module.exclude_sites_with_high_mean(self.args.maxmean)
        if self.args.minstd is not None:
            self.module.remove_lowest_std_sites(self.args.minstd)

        if self.args.rmxy:
            logging.info("excluding sites from X and Y chromosomes...")
            self.module.exclude(loadtxt(HUMAN_X_Y, dtype = str))
        if self.args.rmns:
            logging.info("excluding non-specific sites...")
            self.module.exclude(loadtxt(NONSPECIFIC_PROBES, dtype = str))
        if self.args.rmpoly:
            logging.info("excluding polymorphic sites...")
            self.module.exclude(loadtxt(POLYMORPHIC_CPGS, dtype = str))

    # must  be called after run
    def preprocess_samples_data(self):
        if self.args.keep is not None: # important check, otherwise will keep [] samples (will remove everything)
            self.module.keep(self.keep_list)
        if self.args.remove is not None:
            self.module.remove(self.remove_list)
        if self.args.maxpcstd is not None:
            self.module.exclude_maxpcstds(self.args.maxpcstd)

    # must be called after all preprocessing (preprocess_samples_data, preprocess_sites_data)
    # save methylation data in Glint format
    def save(self, output_perfix):
        if self.args.gsave:
            self.module.save_serialized_data(output_perfix)
        if self.args.txtsave:
            self.module.save_raw_data(output_perfix)

    def run(self, args):
        try:
            self.args = args
            self.module = None
            if args.datafile.name.endswith(methylation_data.GLINT_FILE_SUFFIX):
                logging.info("Loading glint file: %s..." % args.datafile.name)
                self.module = load(args.datafile) # datafile is fileType (status: open for read)
                logging.debug("Got methylation data with %s sites and %s samples id" % (self.module.sites_size, self.module.samples_size))
                # if phenotype or covariates supplied with metylation data, replace module covar and pheno file with new ones
                if args.phenofile is not None:
                    self.module.add_pheno_files(args.phenofile)
                if args.covarfile is not None:
                    self.module.add_covar_files(args.covarfile)
            else:
                self.module = methylation_data.MethylationDataLoader(datafile = args.datafile, phenofile = args.phenofile, covarfiles = args.covarfile)

            # load remove/keep sites/samples files and remove/keep values
            self.include_list = []
            self.exclude_list = []
            self.remove_list = []
            self.keep_list = []

            if args.include is not None:
                self.include_list = self._load_and_validate_ids_in_file(args.include, self.module.cpgnames)
            if args.exclude is not None:
                self.exclude_list = self._load_and_validate_ids_in_file(args.exclude, self.module.cpgnames)
            if args.keep is not None:
                self.keep_list = self._load_and_validate_ids_in_file(args.keep, self.module.samples_ids)
            if args.remove is not None:
                self.remove_list = self._load_and_validate_ids_in_file(args.remove, self.module.samples_ids)

        except Exception:
            logging.exception("in methylation data")
            raise
