import logging
from pickle import load
from modules import methylation_data
from utils import common

GLINT_FORMATTED_EXTENSION = ".glint" #TODO move to a config file

class MethylationDataParser(object):

    ALL_ARGS = ['--datafile', '--include', '--exclude', '--keep', '--remove', '--excludemin', '--excludemax', '--gsave']

    @staticmethod
    def init_args( parser ):
        parser.add_argument('--datafile', required=True, help = "A data matrix file of beta-normalized methylation levels or a .glint file")
    
        group1 = parser.add_mutually_exclusive_group(required = False)
        group1.add_argument('--include', type = str,  help = "A list of sites to include in the data; removes the rest of the sites")
        group1.add_argument('--exclude', type = str,   help = "A list of sites to exclude from the data; includes the rest of the sites")

        group2 = parser.add_mutually_exclusive_group(required = False)
        group2.add_argument('--keep',   type = str, help = "A list of samples to include in the data; removes the rest of the samples")
        group2.add_argument('--remove', type = str, help = "A list of samples to exclude in the data; includes the rest of the samples")

        parser.add_argument('--excludemin', type = float, help = "A threshold for the minimal mean methylation level to consider")
        parser.add_argument('--excludemax', type = float, help = "A threshold for the maximal mean methylation level to consider")

        parser.add_argument('--gsave', action='store_true', help = "Save the data in a glint format; makes following executions faster")

    def _validate_file_exists(self, filepath):
        if not os.path.exists(filepath) :
            logging.error("The file '%s' doesn't exist. Exiting" % filepath)
            common.terminate(self.__class__.__name__)

    def _load_and_validate_file_of_dimentions(self, filepath, dim):
        """
        validates that a file exists and that it is a vector of dimentions dim
        """
        self._validate_file_exists(filepath)
        logging.info("Loading file %s..." % filepath)
        data = loadtxt(filepath, dtype = str)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if len(data.shape) != dim:
            logging.error("The file '%s' is not a %sd vector" % (filepath, dim))
            common.terminate(self.__class__.__name__)

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
            logging.warning("The file %s contains ids that are not found in the datafile: %s" % diff)

        return data

    def _load_sites_file(self, filepath):
        """
        loads file contianing cpgnames list
        warns if there are duplicate sites or sites that are not found in the methylation data
        must be called after self.data is initialized
        """
        return self._load_and_validate_ids_in_file(filepath, self.data.cpgnames)

    
    def _load_sample_ids_file(self, filepath):
        """
        loads file contianing sample ids list
        warns if there are duplicate sample ids or samples ids that are not found in the methylation data
        must be called after self.data is initialized
        """
        return self._load_and_validate_ids_in_file(filepath, self.data.samples_ids)

    def _validate_methylation_value(self, num):
        if not (num <= 1 and num >=0):
            logging.error("must be a standard methylation value (float number between 0 and 1)")

    def _validate_min_and_max_mean_values(self, min_value, max_value):
        """
        validates that the min_value is not greater than max_value
        must be called after self.min_value and self.max_value are set
        """
        if min_value is not None and max_value is not None:
            if max_value <= min_value:
                logging.error("min value %s is greater than max value %s" % (min_value, max_value))
                common.terminate(self.__class__.__name__) 


    def __init__(self, args, output_perfix = ""):
        try:
            # validate arguments
            if args.include is not None:
                self._validate_file_exists(args.include) # TODO seperate file exsits and after module load than validate cpgnames
            if args.exclude is not None:
                self._validate_file_exists(args.exclude)
            if args.keep is not None:
                self._validate_file_exists(args.keep)
            if args.remove is not None:
                self._validate_file_exists(args.remove) 
            if args.excludemin is not None:
                self._validate_methylation_value(args.excludemin)
            if args.excludemax is not None:
                self._validate_methylation_value(args.excludemax) 
            self._validate_min_and_max_mean_values(args.excludemin, args.excludemax)

            # read methylation datafile
            if args.datafile.endswith(GLINT_FORMATTED_EXTENSION):
                with open(args.datafile,'rb') as f:
                    logging.info("Loading glint file: %s..." % args.datafile)
                    self.data = load(f)
                    logging.debug("Got methylation data with %s sites and %s samples id" % (self.module.sites_size, self.module.samples_size))
            else:
                self.data = methylation_data.MethylationData(datafile = args.datafile)

            # load remove/keep sites/samples files and remove/keep values
            if args.include is not None:
                include_list = self._load_sites_file(args.include)
                self.data.include(include_list)
            if args.exclude is not None:
                exclude_list = self._load_sites_file(args.exclude)
                self.data.exclude(exclude_list)
            if args.keep is not None:
                keep_list = self._load_sample_ids_file(args.keep)
                self.data.keep(keep_list)
            if args.remove is not None:
                remove_list = self._load_sample_ids_file(args.remove)
                self.data.remove(remove_list)

            # exclude min/max values
            if args.excludemin is not None:
                self.data.exclude_sites_with_low_mean(args.excludemin)
            if args.excludemax is not None:
                self.data.exclude_sites_with_high_mean(args.excludemax)
            
            # save methylation data in Glint format
            if args.gsave:
                self.data.save(output_perfix + methylation_data.COMPRESSED_FILENAME + GLINT_FORMATTED_EXTENSION)

        except Exception:
            logging.exception("in methylation data")
            raise
