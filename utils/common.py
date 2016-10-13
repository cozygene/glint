from time import time
import sys
import logging
import inspect
from pandas import read_csv, DataFrame
from numpy import float32, array

DELIMITERS = {'\t': 'tab',
              ' ' : 'space', 
              ',' : 'comma'}
              
def terminate(error_msg):
    logging.error("ERROR: " + error_msg)
    sys.exit(2)

def get_dim(vector):
    if vector.ndim == 1 or (vector.ndim == 2 and vector.shape[1] == 1):
        return 1
    return 2

def get_delimiter(line):
    lengths = []
    for sep in DELIMITERS.keys():
        lengths.append(len(line.split(sep)))
    max_len = max(lengths) # the right delimiter will split the line to more than one cell
    assert lengths.count(max_len) == 1 # only one possible delimiter
    max_i = lengths.index(max_len)
    for l in lengths:
        assert l in [max_len, 1]
    return DELIMITERS.keys()[max_i]


def loadtxt(filepath, dtype = None, header = None, delimiter='\t'):
    if dtype:
        x = read_csv(filepath, dtype = dtype, header = header, sep=delimiter)
    else:
        x = read_csv(filepath, dtype = object, header = None, sep=delimiter)
    return DataFrame.as_matrix(x)

def load_float_data_and_headers(filepath, delimiter='\t', dtype = float32, na_values = None):
    """
    assumes data is 2 dimension matrix
    gets data matrix (from type str)
    extracts headers if exist and convert data type to float
    returns (None,None,None) if error occured 

    returns data - data without headers (data type is float)
            cols_title - list of the cols header (first line) (None if no header was found)
            rows_title - list of the rows header (first column) (None if no header was found)
    """
    #check if the first line is an header
    col_names = None
    with open(filepath, 'r') as f:
        first_row = f.readline()
    f.close()

    #find right delimiter
    sep = get_delimiter(first_row)
    if sep != delimiter:
        logging.info("Switching to %s delimited matrix..." % DELIMITERS[sep])
        delimiter = sep
        
    first_row = first_row.strip().split(delimiter)
    num_of_cols = len(first_row)


    try:
        float(first_row[-1])
    except ValueError: # this is a header and not part of the data
        col_names = first_row[1:]

    # check if the first column is a row titles
    row_names = None
    first_col = DataFrame.as_matrix(read_csv(filepath, dtype=str, delimiter=delimiter, usecols=[0], header=None))
    first_col = first_col.reshape(-1,)
    try:
        float(first_col[-1])
    except ValueError: # this is a row's title and not part of the data
        if col_names:
            row_names = first_col[1:]
        else:
            row_names = first_col
    # read the data - if there are row names, ignore first column. otherwise read the whole data    
    if col_names:   # there is an header in the data
        header = 0      # pandas will infer it
    else:           # there is no header in the data
        header = None

    try:
        if row_names is not None:
            if col_names is not None:
                data = read_csv(filepath, dtype=dtype, delimiter=delimiter, header=header, usecols=col_names, na_values=na_values)
            else:
                data = read_csv(filepath, dtype=dtype, delimiter=delimiter, header=header, usecols=range(1,num_of_cols), na_values=na_values)
        else:
            data = read_csv(filepath, dtype=dtype, delimiter=delimiter, header=header, na_values=na_values)
        
        
    except Exception as e:
        logging.exception("While loading data")
        terminate("Error with the format of the file. delimiter must be one of %s and float data" % ", ".join(DELIMITERS.keys())) 

    # convert to array
    if col_names is not None:
        col_names = array(col_names) # by default will be 1d
    if data is not None:
        data = DataFrame.as_matrix(data)
    else:
        common.terminate("No data were found in the file %s." % filepath)

    return data, col_names, row_names

def load_data_file(filepath, dim, dtype = float32, na_values = None):
    """
    filepath is the path to the file to load
    dim is the minimal dimension of the matrix in the file
    """
    data = None
    a = time()
    try:
        data, col_names, row_names = load_float_data_and_headers(filepath, delimiter='\t', dtype = dtype, na_values = na_values)
    except Exception as e:
        logging.exception("While loading data")
        terminate("Some error with the data file format; please check the delimiter.")
    
    logging.debug("Read with pandas took %s seconds." % (time()-a))
    
    if data is None:
        terminate("Could not read data file.")
    if data.size == 1:
        logging.warning("Only one value found in the file: %s." % data)

    return data, col_names, row_names
