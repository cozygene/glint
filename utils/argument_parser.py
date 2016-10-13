import argparse
import sys


class GlintMutuallyExclusiveGroup(argparse._MutuallyExclusiveGroup):
        def __init__(self, container, required=False):
            super(GlintMutuallyExclusiveGroup, self).__init__(container, required)
            self._all_args = set()
            self._required_arguments = []
            self._arguments_dependencies = {}

        def add_mutually_exclusive_group(self, **kwargs):
            raise Exception("NOT SUPPORTED: func add_mutually_exclusive_group in GlintMutuallyExclusiveGroup")

        def argname(self, arg_str):
            return arg_str.replace('-','')

        def get_required_args(self):
            return self._required_arguments

        # TODO check if duplicate flags are identical
        def add_argument( self,  *args, **kwargs ):
            arg = self.argname(args[0])

            if 'required' in kwargs:
                self._required_arguments.append(arg)
                kwargs['required'] = False

            if 'dependencies' in kwargs:
                self._arguments_dependencies[arg] = [self.argname(dependency) for dependency in kwargs.pop('dependencies')]

            self._all_args.add(args[0])
    

            for option_string in args:
                if option_string in self._option_string_actions: #TODO check that they are equal
                    self._option_string_actions.pop(option_string) 
            # this will reenter the last one and wont delete it from description
            super(GlintMutuallyExclusiveGroup, self).add_argument(*args, **kwargs)

        def get_all_args(self):
            return self._all_args

        def get_args_dependencies(self):
            return self._arguments_dependencies


class GlintArgumentGroup(argparse._ArgumentGroup):
        
        def __init__(self, container, title=None, description=None, **kwargs):
            super(GlintArgumentGroup, self).__init__(container, title, description, **kwargs)

            # glint additional arguments per group
            self._all_args = set()
            self._required_arguments = []
            self._arguments_dependencies = {}
            self._glint_mutually_exclusive_groups = [] # must not use _mutually_exclusive_groups since _mutually_exclusive_groups is common for all mutually exclusive groups in all groups


        def add_mutually_exclusive_group(self, **kwargs):
            group = GlintMutuallyExclusiveGroup(self, **kwargs)
            self._mutually_exclusive_groups.append(group)
            self._glint_mutually_exclusive_groups.append(group)
            return group

        def argname(self, arg_str):
            return arg_str.replace('-','')

        # TODO check if duplicate flags are identical
        def add_argument(self, *args, **kwargs):
            arg = self.argname(args[0])

            if 'required' in kwargs:
                self._required_arguments.append(arg)
                kwargs['required'] = False
            
            if 'dependencies' in kwargs:
                self._arguments_dependencies[arg] = [self.argname(dependency) for dependency in kwargs.pop('dependencies')]

            self._all_args.add(args[0])

            for option_string in args:
                if option_string in self._option_string_actions: #TODO check that they are equal
                   self._option_string_actions.pop(option_string) 

            # this will reenter the last one and wont delete it from description
            super(GlintArgumentGroup, self).add_argument(*args, **kwargs)


        def get_required_args(self):
            required_args = self._required_arguments[:] # copy by value, not by reference
            [required_args.extend(mutual_group.get_required_args()) for mutual_group in self._glint_mutually_exclusive_groups]
            return required_args

        def get_args_dependencies(self):
            for group in self._glint_mutually_exclusive_groups :
                self._arguments_dependencies.update(group.get_args_dependencies())
            return self._arguments_dependencies

        def get_all_args(self):
            for group in self._glint_mutually_exclusive_groups :
                self._all_args.update(group.get_all_args())
            return self._all_args

class GlintArgumentParser(argparse.ArgumentParser):#,argparse._ActionsContainer, GlintArgumentGroup):

    def __init__(self, *args, **kwargs):
        super(GlintArgumentParser, self).__init__(*args, **kwargs)
        self.all_modules = []

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

    def add_argument_group(self, *args, **kwargs):
        group = GlintArgumentGroup(self, *args, **kwargs)
        self._action_groups.append(group)
        return group

    def add_mutually_exclusive_group(self, **kwargs):
        raise Exception("NOT SUPPORTED: func add_mutually_exclusive_group in GlintArgumentParser")

    def add_argument( self,  *args, **kwargs ):
        raise Exception("NOT SUPPORTED: func add_argument in GlintArgumentParser")

            