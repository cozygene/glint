import sys
import logging
import inspect

def terminate(module_name):
    logging.info("%s was terminated" % module_name)
    sys.exit(2)