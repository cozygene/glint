import argparse
import os
import sys
from numpy import loadtxt

from parsers import ewas
from parsers import refactor


MODULES = ["ewas", "refactor"] 

class GlintArgumentParser(argparse.ArgumentParser):
    # TODO add epilog
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)


def init_parser( parser ):
    parser.add_argument('--datafile', required=True, help = "path to a data matrix of beta-normalized methylation levels")
    parser.add_argument('--refactor', action='store_true', help = "help") 
    parser.add_argument('--ewas',      action='store_true', help = "help") 

# TODO make it be posiible to run from import
def parse ( args ):
    if args.datafile and not os.path.exists(args.datafile) :
        logging.error("The file '%s' doesn't exist. Exiting" % args.datafile)
        sys.exit(2)
    
    data = loadtxt(args.datafile, dtype = str)
    sample_ids = data[0,:][1:]
    methylation_sites = data[:,0][1:]
    data = data[1:,1:].astype(float)

    if args.refactor:
        refactor.parse(data, args)

    if args.ewas:
        ewas.parse(args)

if __name__ == '__main__':

    parser = GlintArgumentParser(prog=os.path.basename(sys.argv[0]), conflict_handler='resolve', description = "<< add help before >>", epilog = "<< add help after >>")
    
    # set arguments
    init_parser(parser)
    for m in MODULES:
        globals()[m].init_parser( parser )
    
    # parse arguments
    # TODO should we call any module parse?
    args = parser.parse_args() 
    parse(args)




