import argparse
import sys

class GlintMutuallyExclusiveGroup(argparse._MutuallyExclusiveGroup):
        def add_argument( self,  *args, **kwargs ):
            for option_string in args:
                if option_string in self._option_string_actions: #TODO check that they are equal
                    old = self._option_string_actions.pop(option_string) #will remove this from thew list so argparse won;t except on duplicate
                    
            return super(GlintMutuallyExclusiveGroup, self).add_argument( *args, **kwargs)

            

class GlintArgumentGroup(argparse._ArgumentGroup):
        #TODO check if equal 

        # print argparse._StoreAction(**kwargs) == self._option_string_actions[option_string]

        def add_argument( self,  *args, **kwargs ):
            for option_string in args:
                if option_string in self._option_string_actions: #TODO check that they are equal
                    self._option_string_actions.pop(option_string) 
            # doc: this will reenter the last one and wont delete it from description
            return super(GlintArgumentGroup, self).add_argument( *args, **kwargs)


            

class GlintArgumentParser(argparse.ArgumentParser):#,argparse._ActionsContainer, GlintArgumentGroup):

    # TODO add epilog

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

    def add_argument_group(self, *args, **kwargs):

        group = GlintArgumentGroup(self, *args, **kwargs)
        self._action_groups.append(group)
        return group

    def add_mutually_exclusive_group(self, **kwargs):
        print "!@#$!$!@#$!@$!@$!@#$#!$!@$#!$!#!#$#$##$$#"
        group = GlintMutuallyExclusiveGroup(self, **kwargs)
        self._mutually_exclusive_groups.append(group)
        return group

    def add_argument( self,  *args, **kwargs ):
        for option_string in args:
            if option_string in self._option_string_actions: #TODO check that they are equal
                old = self._option_string_actions.pop(option_string) #will remove this from thew list so argparse won;t except on duplicate
                
        return super(GlintArgumentParser, self).add_argument( *args, **kwargs)

            