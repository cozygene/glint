#!/usr/bin/env python
import argparse
import os
import sys
from configuration import configurelogging
configurelogging.configureLogging('') #todo should seperate each module to a different folder to have different "namespaces"?
import logging
from utils import common
from numpy import loadtxt
from utils import GlintArgumentParser
from parsers import ModuleParser, RefactorParser, EWASParser, MethylationDataParser  #dont remove this is imported in,,,

class GlintParser(ModuleParser):
    def __init__(self, parser):
        optional = parser.add_argument_group('3.Optional arguments')
        optional.add_argument('-h', '--help', action='help', help = "print this help") # add help here so it will be under same group with all other optional argument
        optional.add_argument('--out', type = str,  default = "",   help = "changes the prefix of the output file ")

        modules = parser.add_argument_group('4.Glint modules')
        modules.add_argument('--refactor', action='store_true', help = "<todo add help here>")
        modules.add_argument('--ewas',     action='store_true', help = "<todo add help here>" )

        super(GlintParser, self).__init__(optional, modules)
    

class ModulesArgumentParsers(object):
    FUNCTIONALITY_ARGS = [ '--gsave', '--refactor', '--ewas'] # TODO find better way to hold arguments that cause some functionality. glint is not supposed to be aware of those args
    DATA_PREPROCESSING_NOT_RELEVANT_FOR_REFACTOR = ['--include', '--exclude', '--minmean', '--maxmean']
        
    def __init__(self, user_args_selection):
        self.selected_args = user_args_selection
        self.parser = GlintArgumentParser(prog=os.path.basename(sys.argv[0]),
                             description = "<<todo add help before >>",
                             epilog = "<< todo add help after >>",
                             add_help=False) # don't add help because it is added in 'optional group'

        self.glint_parser = None
        self.meth_parser = None
        self.refactor_parser = None
        self.ewas_parser = None
        self.args = None

    def add_arguments(self):
        # Notice to add arguments by the order you want them to be printed in --help

        #main
        self.meth_parser = MethylationDataParser(self.parser)
        self.glint_parser = GlintParser(self.parser)

        #modules
        self.refactor_parser = RefactorParser(self.parser)
        self.ewas_parser = EWASParser(self.parser)

    def parse_args(self):
        logging.info("Validating arguments...")
        self.args = self.parser.parse_args()

        optional_args = []

        self.meth_parser.validate_args(self.args)
        optional_args.extend(self.meth_parser.all_args)

        self.glint_parser.validate_args(self.args)
        optional_args.extend(self.glint_parser.all_args)
        
        if self.args.refactor:
            self.refactor_parser.validate_args(self.args)
            optional_args.extend(self.refactor_parser.all_args)
        if self.args.ewas:
            self.ewas_parser.validate_args(self.args)
            optional_args.extend(self.ewas_parser.all_args)

        self.check_selected_args(optional_args)
        return self.args


    def check_selected_args(self, optional_args):
        """
        validates that the user selected "action" argument (and didn't supply just a datafile) 
        and that the user didnt select an argument that is not an option for him (argument from a module that wasn't selected)
        """
        selected_args = set(self.selected_args)
        func_args = set(self.FUNCTIONALITY_ARGS)
        optional_args = set(optional_args)
        args_not_relevant_for_refactor = set(self.DATA_PREPROCESSING_NOT_RELEVANT_FOR_REFACTOR)

        if len(selected_args.intersection(func_args)) == 0:
            common.terminate("Nothing to do with the data, select at least one argument from %s" % self.FUNCTIONALITY_ARGS)

        differ = selected_args.difference(optional_args)
        if differ:
            common.terminate("selected redundent argument" + str(differ)) # TODO: tell which module the arguments belong to, to user might forgot to specify --module and that msg can be confusing

        # warn if only refactor module is selected and  user selectes data management flags that are not relevant for refactor
        if self.args.refactor and len(selected_args.intersection(func_args)) == 1:
            not_relevant_atgs = list(selected_args.intersection(args_not_relevant_for_refactor))
            if len(not_relevant_atgs) != 0:
                logging.warning("selected data management arguments which are not relevant for refactor: %s" % str(not_relevant_atgs))

    def run(self):
        self.meth_parser.run(self.args)
        self.meth_parser.preprocess_samples_data() # preprocess samples before refactor and before ewas

        if self.args.refactor:
            refactor_meth_data = self.meth_parser.module.copy()
            self.refactor_parser.run(args = self.args,
                                    meth_data = refactor_meth_data,
                                    output_perfix = self.args.out)
            self.meth_parser.module.add_covariates(self.refactor_parser.module.components) # add refactor components as covariate file

        self.meth_parser.preprocess_sites_data() #preprocess sites after refactor and before ewas
        self.meth_parser.gsave(output_perfix = self.args.out) #save after all preprocessing #TODO maybe take gsave out from MethData module

        if self.args.ewas:
            ewas_meth_data = self.meth_parser.module.copy()
            self.ewas_parser.run(args = self.args,
                                 meth_data = ewas_meth_data,
                                 output_perfix = self.args.out)

if __name__ == '__main__':
    selected_args = [arg for arg in sys.argv if arg.startswith("-")] #TODO startwith "--"? (there are no arguments that starts with -)
    
    parser = ModulesArgumentParsers(selected_args)

    parser.add_arguments()
    args = parser.parse_args()

    logging.info("Starting glint...")
    parser.run()




