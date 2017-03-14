import logging
from utils import common

"""
HOW TO write a new parser:
1. create new class that inherits from ModuleParser
2. in the __init__ function define all the arguments of this module, note that:
    a. you first should declare on new argument group (or groups) (with add_argument_group). you can see the name you give when using --help
    b. declare all arguments. associate them with the group 
    c. You can declare required arguments ONLY if you associated the arguments with a uniqe group for this module (otherwise glint will require them any time and not only when module is selected)
    d. arguments must begin with "--" (and not "-")
    e. you can define dependencies list for each argument with the flag 'dependencies' (e.g dependencies = ["--pheno"]).
      This will validate that if the user choose this argument, he must also choose the arguments in the dependencies list (e.g '--pheno')
3. if you need to validate condition on your argument before passing them to the real module to handle them - do that in function validate_args
4. override function 'run' which should call and run your module there

In case when your parser depend on another parser (e.g PlotParser depend on QQPlotParser and ManhattanPlotParser)
    1. initiates the parsers in the __init__ function and save them as members
        (self.dependency_parser = DependencyParser(parser) )
    2. add the parser flags and required flags to your parser flags:
            self.all_args.extend(self.dependency_parser.all_args)
            self.required_args.extend(self.dependency_parser.required_args)

    * if you dont always "need your dependencies" (i.e qq-plot cant run with manhattan plot)
        you can add to your parser a flag for each of it's  dependencies (i.e PlotParser has a --qqplot flag  and --manhattan flag)
        and on your validate_args function - initiate the needed parser if the right flag was selected.
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
        self.required_args = []
        for group in self.groups:
            self.all_args.extend(group.get_all_args())
            self.required_args.extend(group.get_required_args())


    def validate_required_args(self, args):
        """
        errors if required argument is not set in arguments list (args)
        """
        for arg in self.required_args:
            if not contains_arg(args, arg):
                common.terminate("argument --%s is required" % arg)

    def validate_args_dependencies(self, args):
        for group in self.groups:
            self._validate_args_dependencies_in_group(group, args)


    def _validate_args_dependencies_in_group(self, group, args):
        """
        errors if there is an argument in argument list (args) that is set but its dependent argument are not set int that list
        """
        for arg, dependencies in group.get_args_dependencies().iteritems():
            if contains_arg(args, arg):
                for dependency in dependencies:
                    if not contains_arg(args, dependency):
                        common.terminate("argument --%s requires argument --%s" % (arg, dependency))


    def validate_args(self, args):
        self.validate_required_args(args)
        self.validate_args_dependencies(args)
    
