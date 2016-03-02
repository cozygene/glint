import logging
from utils import common


def contains_arg(args, arg):
    """
    checks if argument (arg) is set in arguments list (args)
    """
    return args.__getattribute__(arg) is not None

class ModuleParser(object):

    def __init__(self,  *groups):
        self.groups = groups # groups of type GlintArgumentGroup
        self.all_args = []
        for group in self.groups:
            self.all_args.extend(group.get_all_args())


    def validate_required_args(self, args):
        for group in self.groups:
            self._validate_required_args_in_group(group, args)


    def _validate_required_args_in_group(self, group, args):
        """
        errors if required argument is not set in arguments list (args)
        """
        for arg in group.get_required_args():
            if not contains_arg(args, arg):
                logging.error("argument --%s is required" % arg)
                common.terminate(self.__class__.__name__)

    def validate_args_dependencies(self, args):
        for group in self.groups:
            self._validate_args_dependencies_in_group(group, args)


    def _validate_args_dependencies_in_group(self, group, args):
        """
        errors if there is an argument in argument list (args) that is set but it's dependent argument are not set int that list
        """
        for arg, dependencies in self.group.get_args_dependencies().iteritems():
            if contains_arg(args, arg):
                for dependency in dependencies:
                    if not contains_arg(args, dependency):
                        logging.error("argument --%s requires argument --%s" % (arg, dependency))
                        common.terminate(self.__class__.__name__)


    def validate_args(self, args):
        self.validate_required_args(args)
        self.validate_args_dependencies(args)
    
