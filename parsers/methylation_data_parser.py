from modules import methylation_data

class MethylationDataParser(object):

    ALL_ARGS = ['--datafile', '--include', '--exclude', '--keep', '--remove']

    @staticmethod
    def init_args( parser ):
        parser.add_argument('--datafile', required=True,      help = "path to a data matrix of beta-normalized methylation levels")
    
        group1 = parser.add_mutually_exclusive_group(required = False)
        group1.add_argument('--include', type = str,  help = "file with a list of cpgs to include in the data (remove the rest of the cpgs)")
        group1.add_argument('--exclude', type = str,   help = "file with a list of cpgs to exclude from the data (keep the rest of the cpgs)")

        group2 = parser.add_mutually_exclusive_group(required = False)
        group2.add_argument('--keep',   type = str, help = "file with a list of samples to include in the data (remove the rest of the samples)")
        group2.add_argument('--remove', type = str, help = "file with a list of samples to exclude from the data (keep the rest of the samples)")

    def __init__(self, args):
        # init module to validate other args format
        self.module  = methylation_data.MethylationData(datafile = args.datafile,
                                                        includefile = args.include,
                                                        excludefile = args.exclude,
                                                        keepfile = args.keep,
                                                        removefile = args.remove)
