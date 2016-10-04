from time import time
import sys
import logging
import inspect
from pandas import read_csv, DataFrame


def terminate(error_msg):
    logging.error("ERROR: " + error_msg)
    sys.exit(2)


def loadtxt(filepath, dtype = float, header = None, delimiter='\t'):
    x = read_csv(filepath, dtype = dtype, header = header, sep=delimiter)
    return DataFrame.as_matrix(x)

def find_and_extract_headers_in_data(data, dim):
    """
    assumes data is 2 dimension matrix
    gets data matrix (from type str)
    extracts headers if exist and convert data type to float
    returns (None,None,None) if error occured 

    returns data - data without headers (data type is float)
            cols_title - list of the cols header (first line) (None if no header was found)
            rows_title - list of the rows header (first column) (None if no header was found)
    """
    cols_title = None
    rows_title = None

    #search for an header (title for each column)
    try:
        data[0,:][1:].astype(float)
    except:
        cols_title = data[0,:] # header (e.g samples ID)

    # search for title for each row
    try:
        data[:,0][1:].astype(float)
    except:
        rows_title = data[:,0] # first column (e.g sites)

    # remove the header from data
    if cols_title is not None and rows_title is not None:
        data = data[1:,1:]
        # remove the cell at (0,0) (the "title of the entire matrix")
        cols_title = cols_title[1:]
        rows_title = rows_title[1:]
    elif cols_title is not None and rows_title is None:
        data = data[1:,:]
    elif cols_title is None and rows_title is not None:
        data = data[:,1:]
    # else (both are None) dont change the data

    # try to convert data to float
    try:
        data = data.astype(float) 
    except ValueError:
        terminate("file contains values which are not float") # todo change when missing values are supported
    
    if rows_title is not None:
        rows_title = rows_title.astype(str)
    
    if cols_title is not None:
        cols_title = cols_title.astype(str)

    return data, cols_title, rows_title

def load_data_file(filepath, dim):
    """
    filepath is the path to the file to load
    dim is the dimension of the matrix in the file
    """
    data = None
    a = time()
    try:

        data = loadtxt(filepath, dtype = str, header = None, delimiter='\t')
    except:
        terminate("some error with the data file format")
    logging.debug("read with pandas took %s seconds" % (time()-a))
    
    if data is None:
        terminate("some error with the data file")
    if data.size == 1:
        logging.warning("only one value found in the file: %s" % data)
    elif data.ndim != dim:
        return None, None, None

    #headers are assumed to be in data of 2 dimensions (and not in a vector of 1D)    
    if dim == 2:
        return find_and_extract_headers_in_data(data, dim)
    # just change data type
    else:
        try:
            data = data.astype(float)
        except ValueError:
            terminate("file contains values which are not float") # todo change when missing values are supported

    return data, None, None
