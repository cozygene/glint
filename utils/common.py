import sys
import logging
import inspect
from numpy import loadtxt

DELIMITERS = [',',':']

def terminate(error_msg):
    logging.error("ERROR: " + error_msg)
    sys.exit(2)

def load_data_file(filepath, dim):
    """
    filepath is the path to the file to load
    dim is the dimension of the matrix in the file
    """    
    data = loadtxt(filepath, dtype = str)
    sep_i = 0
    while data.ndim != dim and sep_i < len(DELIMITERS):
        data = loadtxt(filepath, dtype = str, delimiter=DELIMITERS[sep_i]) # try different line seperator
        sep_i += 1
    
    if data.ndim != dim:
        return None
    
    return data