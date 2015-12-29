import argparse
import os
import sys
import logging
import copy
from numpy import genfromtxt ,loadtxt
from utils import GlintArgumentParser
from parsers import RefactorParser, EWASParser, MethylationDataParser  #dont remove this is imported in,,,


MODULES_PARSERS = ['MethylationDataParser', 'RefactorParser', 'EWASParser'] 

ALL_ARGS = ['--refactor', '--ewas', '--out']
def add_arguments(parser):       
    parser.add_argument('--refactor', action='store_true', help = "help")
    parser.add_argument('--ewas',     action='store_true', help = "help" )
    parser.add_argument('--out',      type = str,  default ="",    help = "changes the prefix of the output file ")

    
    for m in MODULES_PARSERS:
        globals()[m].init_args(parser)


# TODO make it be posiible to run from import
def run ( args ):
    output_file_prefix = args.out
    optional_args = ALL_ARGS
    modules_to_run = []

    optional_args.extend(MethylationDataParser.ALL_ARGS)
    meth_data = MethylationDataParser(args).module

    # init modules (the init verifies the arguments)
    if args.refactor:
        refactor_meth_data = copy.deepcopy(meth_data)
        optional_args.extend(RefactorParser.ALL_ARGS)
        modules_to_run.append(RefactorParser(args, refactor_meth_data))

    if args.ewas:
        optional_args.extend(EWASParser.ALL_ARGS)
        modules_to_run.append(EWASParser(args))
    
    # run every module
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

    differ = set(selected_args).difference(set(optional_args))
    if differ:
        print "ERROR: selected unused argument" + str(differ)




