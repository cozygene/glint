from time import time
import sys
import logging
import inspect
from pandas import read_csv, DataFrame
from numpy import float32, array

DELIMITERS = {'\t': 'tab',
              ' ' : 'space', 
              ',' : 'comma'}

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def substr(s):
    """
    given string S of length n, returns the string S(0)...S(i)
    where S(i+1) is the first char which is not a string  (an number)
    """
    substr=""
    for i in s:
        if not is_number(i):
            substr+=i
        else:
            break
    return substr

def is_there_a_header(first_col):
    """
    we assume the first column of the data file contains the rows labels ("titles")
    therefor, if there is an "ID" cell (at the corner, the first cell of the first column)
    we can conclude that there is an header to the file (first row, column labels).
    so we need to figure if the first "cell" is the first row's label or an "ID"
    """
    first = first_col[0]
    second = first_col[1]

    if is_number(first) != is_number(second):
        # if they are from different types (one is a number and the other is not) - there is an ID
        return True

    elif substr(first) != substr(second):
        # if the begining of each "cell" is different, there is an "ID"
        return True

    return False

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

def has_header(first_line):
    try:
        float(first_line[-1])
    except ValueError: #first cell is string
        return True
    return False

def loadtxt(filepath, dtype = None, header = None, delimiter='\t'):
    if dtype:
        x = read_csv(filepath, dtype = dtype, header = header, sep=delimiter)
    else:
        x = read_csv(filepath, dtype = object, header = None, sep=delimiter)
    return DataFrame.as_matrix(x)

def load_float_data_and_headers(filepath, header=None, delimiter='\t', dtype = float32, na_values = None):
    """
    dim -  is the minimal dimension of the matrix in the file
    header - True - if file contains a header (first row, which values are the column names)
             False - if file doesn't contain an header
             None - let the function infer if there is an header

    Note that rows titles are auumes to exist (rows titles are the first colum, ids for each row)

    assumes data is 2 dimension matrix
    gets data matrix (from type str)
    extracts headers if exist and convert data type to float
    returns (None,None,None) if error occured 

    returns data - data without headers (data type is float)
            cols_title - list of the cols header (first line) (None if no header was found)
            rows_title - list of the rows header (first column) (None if no header was found)
    """
    #check if the first line is an header
    with open(filepath, 'r') as f:
        first_row = f.readline()
        second_row = f.readline()
    f.close()

        # first row is an header
        
    #find right delimiter
    sep = get_delimiter(second_row) #dont use first row cause it could have only one value(not delimiter)
    if sep != delimiter:
        logging.info("Switching to %s delimited matrix..." % DELIMITERS[sep])
        delimiter = sep
        
    first_row = first_row.strip().split(delimiter)
    second_row = second_row.strip().split(delimiter)
    first_col = DataFrame.as_matrix(read_csv(filepath, dtype=str, delimiter=delimiter, usecols=[0], header=None))
    first_col = first_col.reshape(-1,)
    
    if len(first_row) != len(second_row):
        """
        one option:

        samp1 samp2
        cp1 1 2
        cp2 3 4
        """

        if header == False:
            # files from this format must contain an header
            terminate("It is assumed that data file doesn't contain an header but it does")
        
        # No "ID" at the corner
        # we can be sure we have headers to the rows and cols
        col_names = first_row 
        row_names = first_col[1:]
        data = read_csv(filepath, dtype=dtype, delimiter=delimiter, header=0, usecols=range(1,len(first_row)+1), na_values=na_values)

    else:
        """
        2 options (for 2 samples and 2 sites)

        #1
        ID samp1 samp2
        cp1 1 2 
        cp2 3 4 

        #2
        cp1 1  2 
        cp2 3  4

        """
        try:
            # if there is an header
            if header or (header == None and is_there_a_header(first_col)): #first option, has the "ID cell"
                col_names = first_row[1:]
                row_names = first_col[1:]
                data = read_csv(filepath, dtype=dtype, delimiter=delimiter, header=0, usecols=col_names, na_values=na_values)
            # if there is no header
            else:
                col_names = None
                row_names = first_col
                data = read_csv(filepath, dtype=dtype, delimiter=delimiter, header=None, usecols=range(1,len(first_row)), na_values=na_values)
            
        except Exception as e:
            logging.exception("While loading data")
            terminate("Error with the format of the file. delimiter must be one of %s and float data" % ", ".join(DELIMITERS.keys())) 

    # convert to array
    if col_names is not None:
        col_names = array(col_names) # by default will be 1d
    if data is not None:
        data = DataFrame.as_matrix(data)
    else:
        terminate("No data found in the file %s." % filepath)
    
    return data, col_names, row_names



def load_data_file(filepath, dim, header=None, dtype = float32, na_values = None):
    """
    filepath - is the path to the file to load
    dim -  is the minimal dimension of the matrix in the file
    header - True if file contains a header (first row, which values are the column names)
             False if file doesn't contain an header
             None - let the function infer if there is an header

    Note that rows titles are auumes to exist (rows titles are the first colum, ids for each row)
    """
    data = None
    a = time()
    try:
        data, col_names, row_names = load_float_data_and_headers(filepath, header, delimiter='\t', dtype = dtype, na_values = na_values)
    except Exception as e:
        logging.exception("While loading data")
        terminate("Some error with the data file format; please check the delimiter.")
    
    logging.debug("Read with pandas took %s seconds." % (time()-a))
    
    if data is None:
        terminate("Could not read data file.")
    if data.size == 1:
        logging.warning("Only one value found in the file: %s." % data)

    return data, col_names, row_names