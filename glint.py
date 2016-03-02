import argparse
import os
import sys
from configuration import configurelogging
configurelogging.configureLogging('') #todo should seperate each module to a different folder to have different "namespaces"?
import logging
from numpy import genfromtxt ,loadtxt
from utils import GlintArgumentParser
from parsers import RefactorParser, EWASParser, MethylationDataParser  #dont remove this is imported in,,,

MODULES_PARSERS = ['MethylationDataParser', 'RefactorParser', 'EWASParser'] 
ALL_ARGS = ['--refactor', '--ewas', '--out']


def add_arguments(parser): 
    # Notice to add arguments by the order you want them to be printed in --help

    MethylationDataParser.init_args(parser)

    optional = parser.add_argument_group('3.Optional arguments')
    optional.add_argument('-h', '--help', action='help', help = "print this help") # add help here so it will be under same group with all other optional argument
    optional.add_argument('--out', type = str,  default = "",   help = "changes the prefix of the output file ")

    modules = parser.add_argument_group('4.Glint modules')
    modules.add_argument('--refactor', action='store_true', help = "<todo add help here>")
    modules.add_argument('--ewas',     action='store_true', help = "<todo add help here>" )
    
    
    for m in MODULES_PARSERS:
        globals()[m].init_args(parser)

# TODO make it be posiible to run from import
def run ( args, selected_args):
    output_file_prefix = args.out
    optional_args = ALL_ARGS
    modules_to_run = []

    logging.info("Validating methylation data...")
    optional_args.extend(MethylationDataParser.ALL_ARGS)
    meth_data = MethylationDataParser(args).data

    # init modules (the init verifies the arguments)
    if args.refactor:
        logging.info("Validating refactor...")
        refactor_meth_data = meth_data.copy()
        optional_args.extend(RefactorParser.ALL_ARGS)
        modules_to_run.append(RefactorParser(args, refactor_meth_data))

    if args.ewas:
        logging.info("Validating ewas...")
        optional_args.extend(EWASParser.ALL_ARGS)
        modules_to_run.append(EWASParser(args))

    # validate that the user didnt select an argument that is not an option
    differ = set(selected_args).difference(set(optional_args))
    if differ:
        logging.error("selected redundent argument" + str(differ))

    # run every module
    for m in modules_to_run:
        m.run()


if __name__ == '__main__':
    parser = GlintArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 description = "<<todo add help before >>",
                                 epilog = "<< todo add help after >>",
                                 add_help=False) # don't add help because it is added in 'optional group'
    
    add_arguments(parser)

    # parse arguments
    # TODO should we call any module parse?
    selected_args = [arg for arg in sys.argv if arg.startswith("-")] #TODO startwith "--"? (there are no arguments that starts with -)

    args = parser.parse_args() # modules validation is here
    run(args, selected_args)





