import argparse
import os
import sys
from numpy import loadtxt
from utils import GlintArgumentParser
from parsers import RefactorParser, EWASParser  #dont remove this is imported in,,,


MODULES_PARSERS = ['RefactorParser', 'EWASParser'] 

ALL_ARGS = ['--datafile', '--refactor', '--ewas', '--out']
def add_arguments(parser):       
    parser.add_argument('--datafile', required=True, help = "path to a data matrix of beta-normalized methylation levels")
    parser.add_argument('--refactor', action='store_true', help = "help")
    parser.add_argument('--ewas',      action='store_true', help = "help" )
    parser.add_argument('--out',      type = str,                 help = "changes the prefix of the output file ")

    
    for m in MODULES_PARSERS:
        globals()[m].init_args(parser)


# TODO make it be posiible to run from import
def run ( args ):

    if args.datafile and not os.path.exists(args.datafile) :
        logging.error("The file '%s' doesn't exist. Exiting" % args.datafile)
        sys.exit(2)

    data = loadtxt(args.datafile, dtype = str)
   
    # sample_ids = data[0,:][1:]
    # methylation_sites = data[:,0][1:]
    # data = data[1:,1:].astype(float)

    optional_args = ALL_ARGS
    modules_to_run = []

    if args.refactor:
        optional_args.extend(RefactorParser.ALL_ARGS)
        modules_to_run.append(RefactorParser(args, data))

    if args.ewas:
        optional_args.extend(EWASParser.ALL_ARGS)
        modules_to_run.append(EWASParser(args))
    
    for m in modules_to_run:
        m.run()

    return optional_args


if __name__ == '__main__':

    parser = GlintArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 description = "<< add help before >>",
                                 epilog = "<< add help after >>")# conflict_handler='resolve')
    
    add_arguments(parser)

    # parse arguments
    # TODO should we call any module parse?
    selected_args = [arg for arg in sys.argv if arg.startswith("-")]

    args = parser.parse_args() # modules validation is here

    optional_args = run(args)

    print "-----------"
    print set(selected_args)
    print set(optional_args)
    differ = set(selected_args).difference(set(optional_args))
    if differ:
        print "ERROR: selected unused argument" + str(differ)




