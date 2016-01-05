import logging
from modules import methylation_data

class MethylationDataParser(object):

    ALL_ARGS = ['--datafile', '--include', '--exclude', '--keep', '--remove', '--excludemin', '--excludemax', '--gsave']

    @staticmethod
    def init_args( parser ):
        parser.add_argument('--datafile', required=True,      help = "path to a data matrix of beta-normalized methylation levels or to the Glint file")
    
        group1 = parser.add_mutually_exclusive_group(required = False)
        group1.add_argument('--include', type = str,  help = "file with a list of cpgs to include in the data (remove the rest of the cpgs)")
        group1.add_argument('--exclude', type = str,   help = "file with a list of cpgs to exclude from the data (keep the rest of the cpgs)")

        group2 = parser.add_mutually_exclusive_group(required = False)
        group2.add_argument('--keep',   type = str, help = "file with a list of samples to include in the data (remove the rest of the samples)")
        group2.add_argument('--remove', type = str, help = "file with a list of samples to exclude from the data (keep the rest of the samples)")

        parser.add_argument('--excludemin', type = float, help = "filter out sites with mean methylation level lower than this number. float number between 0 and 1")
        parser.add_argument('--excludemax', type = float, help = "filter out sites with mean methylation level greater than this number. float number between 0 and 1")

        parser.add_argument('--gsave', action='store_true', help = "set to save the methylation data file in Glint format. makes next run faster")

    def __init__(self, args, output_perfix = ""):
        # init module to validate other args format
        try:
            if args.datafile.endswith(".glint"):
                import pickle
                with open(args.datafile,'rb') as f:
                    self.module  = pickle.load(f)
            else :
                self.module  = methylation_data.MethylationData(datafile = args.datafile,
                                                        includefile = args.include,
                                                        excludefile = args.exclude,
                                                        keepfile = args.keep,
                                                        removefile = args.remove,
                                                        min_value = args.excludemin,
                                                        max_value = args.excludemax,
                                                        glint_data_filename = output_perfix + methylation_data.COMPRESSED_FILENAME if args.gsave else None  )
        except Exception:
            logging.exception("in methylation data")
            raise