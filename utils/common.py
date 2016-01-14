import sys
import logging
import inspect

def terminate(module_name):
    import pdb
    pdb.set_trace()
    logging.info("%s was terminated" % module_name)
    sys.exit(2)