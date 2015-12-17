import argparse
import sys
class GlintMutuallyExclusiveGroup(argparse._MutuallyExclusiveGroup):
        def add_argument( self,  *args, **kwargs ):
            print "!@#$!$!@#$!@$!@$!@#$#!$!@$#!$!#!#$#$##$$#"
            for option_string in args:
                if option_string in self._option_string_actions: #TODO check that they are equal
                    old = self._option_string_actions.pop(option_string) #will remove this from thew list so argparse won;t except on duplicate
                    
            return super(GlintMutuallyExclusiveGroup, self).add_argument( *args, **kwargs)

            

class GlintArgumentGroup(argparse._ArgumentGroup):
        #TODO check if equal 
        # print type(self._option_string_actions[option_string])
        # print self._option_string_actions[option_string]
        # print "____"
        # print kwargs
        # print argparse._StoreAction(**kwargs) == self._option_string_actions[option_string]
        # print "HELLO!!!!2"

        # self.register('action', None, _StoreAction)
        # self.register('action', 'store', _StoreAction)
        # self.register('action', 'store_const', _StoreConstAction)
        # self.register('action', 'store_true', _StoreTrueAction)
        # self.register('action', 'store_false', _StoreFalseAction)
        # self.register('action', 'append', _AppendAction)
        # self.register('action', 'append_const', _AppendConstAction)
        # self.register('action', 'count', _CountAction)
        # self.register('action', 'help', _HelpAction)
        # self.register('action', 'version', _VersionAction)
        def add_argument( self,  *args, **kwargs ):
            res = None
            conf = None
            for option_string in args:
                if option_string in self._option_string_actions: #TODO check that they are equal

                    conf = self._option_string_actions[option_string]
          
                    old = self._option_string_actions.pop(option_string) #reutreutreut#will remove this from thew list so argparse won;t except on duplicate
   
            try:    
                res = super(GlintArgumentGroup, self).add_argument( *args, **kwargs)
            except:
                pass

            # print self._group_actions
            # if conf:
                # self._group_actions.append(conf) reutreutreut

                # print self._group_actions
            return res

            

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

            