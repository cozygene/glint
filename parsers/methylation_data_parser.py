import logging
from modules import methylation_data

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