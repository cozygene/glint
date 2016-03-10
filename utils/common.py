import sys
import logging
import inspect

def terminate(error_msg):
    logging.error(error_msg)
    sys.exit(2)