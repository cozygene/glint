import sys
import logging
import inspect

def terminate(error_msg):
    logging.error("ERROR: " + error_msg)
    sys.exit(2)