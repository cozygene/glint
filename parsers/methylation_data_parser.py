import logging
from pickle import load
from modules import methylation_data
from utils import common
import argparse
from module_parser import ModuleParser
GLINT_FORMATTED_EXTENSION = ".glint" #TODO move to a config file

class MethylationDataParser(ModuleParser): 

    def __init__(self, parser):
        required = parser.add_argument_group('1.Required arguments') # numbering in the group name because help print it by abc order
        required.add_argument('--datafile', type = argparse.FileType('r'), required = True,  help = "A data matrix file of beta-normalized methylation levels or a .glint file")

        optional = parser.add_argument_group('2.Data management options ')
        optional.add_argument('--covar',   type = argparse.FileType('r'), nargs='*', help = "A covariates file")
        optional.add_argument('--pheno',   type = argparse.FileType('r'), help = "A phenotype file")
        optional.add_argument('--maxpcstd', metavar=('PC_INDEX (TODO Elior, change those names?', 'STD_COUNT'), type = int, action = 'append', nargs = 2, help = "TODO Elior, edit: pc index and std number of times for removing outliers")
            
        group1 = optional.add_mutually_exclusive_group(required = False)
        group1.add_argument('--include', type = argparse.FileType('r'), help = "A list of sites to include in the data; removes the rest of the sites")
        group1.add_argument('--exclude', type = argparse.FileType('r'), help = "A list of sites to exclude from the data; includes the rest of the sites")

        group2 = optional.add_mutually_exclusive_group(required = False)
        group2.add_argument('--keep',   type = argparse.FileType('r'), help = "A list of samples to include in the data; removes the rest of the samples")
        group2.add_argument('--remove', type = argparse.FileType('r'), help = "A list of samples to exclude in the data; includes the rest of the samples")

        optional.add_argument('--gsave', action='store_true', help = "Save the data in a glint format; makes following executions faster")

        def methylation_value(num):
            num = float(num)
            if not (num <= 1 and num >=0):
                common.terminate("minmean/maxmean must be a standard methylation value (float number between 0 and 1)")
            return num

        optional.add_argument('--minmean', type = methylation_value, help = "A threshold for the minimal mean methylation level to consider")
        optional.add_argument('--maxmean', type = methylation_value, help = "A threshold for the maximal mean methylation level to consider")

        super(MethylationDataParser, self).__init__(required, optional)
        
    def validate_args(self, args):
        super(MethylationDataParser, self).validate_args(args)
        self._validate_min_and_max_mean_values(args.minmean, args.maxmean)

    def _load_and_validate_file_of_dimentions(self, fileobj, dim):
        """
        validates that the file contains a vector of dimentions dim
        """
        if not isinstance(fileobj, file):
            fileobj = open(fileobj, 'r')
        logging.info("loading file %s..." % fileobj.name)
        data = common.load_data_file(fileobj.name, dim)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if data is None:
            common.terminate("The file '%s' is not a %sd vector" % (fileobj.name, dim))

        return data

    def _load_and_validate_ids_in_file(self, filepath, optional_ids_list):
        """
        loads a vector file contianing ids list
        warns if there are duplicate ids in the file or if there are ids which are not found in optional_ids_list
        fails if this file doesn't contains a vector
        """
        data = self._load_and_validate_file_of_dimentions(filepath, 1)
        data_set = set(data)
        if len(data) != len(data_set):
            logging.warning("The file %s contains ids more than once" % filepath)

        diff =  data_set.difference(set(optional_ids_list))
        if diff != set([]):
            logging.warning("The file %s contains ids that are not found in the datafile: %s" % (filepath, diff))

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
    def gsave(self, output_perfix):
        if self.args.gsave:
            if output_perfix:
                output_file_name = output_perfix
            else:
                output_file_name = methylation_data.COMPRESSED_FILENAME
            self.module.save(output_file_name + GLINT_FORMATTED_EXTENSION)

    def run(self, args):
        try:
            self.args = args
            self.module = None
            if args.datafile.name.endswith(GLINT_FORMATTED_EXTENSION):
                logging.info("Loading glint file: %s..." % args.datafile.name)
                self.module = load(args.datafile) # datafile is fileType (status: open for read)
                logging.debug("Got methylation data with %s sites and %s samples id" % (self.module.sites_size, self.module.samples_size))
                # if phenotype or covariates supplied with metylation data, replace module covar and pheno file with new ones
                if args.pheno is not None:
                    self.module.upload_new_phenotype_file(args.pheno)
                if args.covar is not None:
                    self.module.upload_new_covaritates_files(args.covar)
            else:
                self.module = methylation_data.MethylationData(datafile = args.datafile, phenofile = args.pheno, covarfiles = args.covar)

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
