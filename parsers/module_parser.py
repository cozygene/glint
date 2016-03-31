import logging
from utils import common

"""
HOW TO write a new parser:
1. create new class that inherits from ModuleParser
2. in the __init__ function define all the arguments of this mudole, note that:
    a. you first should declare on new argument group (or groups) (with add_argument_group). you can see the name you give when using --help
    b. declare all arguments. associate them with the group 
    c. You can declare required arguments ONLY if you associated the arguments with a uniqe group for this mudole (otherwise glint will require them any time and not only when module is selected)
    d. arguments must begin with "--" (and not "-")
    e. you can define dependencies list for each argument with the flag 'dependencies' (e.g dependencies = ["--pheno"]).
      This will validate that if the user choose this argument, he must also choose the arguments in the dependencies list (e.g '--pheno')
3. if you need to validate condition on your argument before passing them to the real module to handle them - do that on validate_args
4. override function 'run', call and run your module there

TODO - move this to instructions file
"""

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
                common.terminate("argument --%s is required" % arg)

    def validate_args_dependencies(self, args):
        for group in self.groups:
            self._validate_args_dependencies_in_group(group, args)


    def _validate_args_dependencies_in_group(self, group, args):
        """
        errors if there is an argument in argument list (args) that is set but it's dependent argument are not set int that list
        """
        for arg, dependencies in group.get_args_dependencies().iteritems():
            if contains_arg(args, arg):
                for dependency in dependencies:
                    if not contains_arg(args, dependency):
                        common.terminate("argument --%s requires argument --%s" % (arg, dependency))


    def validate_args(self, args):
        self.validate_required_args(args)
        self.validate_args_dependencies(args)
    
