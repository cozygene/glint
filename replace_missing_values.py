from numpy import savetxt, where, zeros, delete, mean, vstack, column_stack
import argparse
from utils import common
from numpy import float32, isnan, nanmean
from numpy import sum as npsum
DATA_TYPE = float32
STR_DATA_TYPE = '|S16'
def _replace_missing_values_in_matrix(all_data, missing_value_indicator, data_max_missing_values, samples_max_missing_values, replace = False, col_names=None, row_names=None):
    number_of_data, number_of_samples = all_data.shape
    na_count_per_sample  = zeros(number_of_samples)
    data_indices_to_remove = []
    print "Replacing missing values by mean..."
    
    na_count_per_sample = npsum(isnan(all_data), axis=0)
    samples_indices_to_keep = where(na_count_per_sample <= samples_max_missing_values*number_of_samples)[0]
    print("%s samples were not replaced because they have more than %s missing values" % (number_of_samples - len(samples_indices_to_keep), samples_max_missing_values))

    na_count_per_site = npsum(isnan(all_data), axis=1)
    sites_indices_to_keep = where(na_count_per_site <= data_max_missing_values*number_of_data)[0]
    print("%s sites were not replaced because they have more than %s missing values" % (number_of_data - len(sites_indices_to_keep), data_max_missing_values))

    print("replacing each missing value by it's site mean...")
    all_data =  all_data[sites_indices_to_keep, :][:,samples_indices_to_keep]
    sites_mean = nanmean(all_data, axis=1)
    for i, data_for_all_samples in enumerate(all_data):
        na_indices = where(isnan(all_data[i,:]))[0]
        all_data[i][na_indices] = sites_mean[i]
   
    # return the relevant samples ids and sites ids
    if col_names is not None:
        col_names = col_names[samples_indices_to_keep]

    if row_names is not None:
        row_names = row_names[sites_indices_to_keep]

    return all_data, col_names, row_names

def get_data_type(data_type_str):
    if data_type_str in ['float', 'double', float]:
        return float
    if data_type_str in ['int', 'integer', int]:
        return int
    return None

def replace_missing(data_filename, missing_value_indicator, data_max_missing_values, samples_max_missing_values, sep = " ", suffix = ".no_missing_values", header = None):
    """
    replaces missing values by mean (mean of non-missing data/samples) and saves the output to the file named data_filename + suffix
    if there are too many missing samples (more than samples_max_missing_values) - they are removed
    if there are too mant missing datas (more than data_max_missing_values) - they are removed

    parames:
    data_filename - a matrix of type int or float (not including the missing value indicator which can be string as well). dimensions nXm where n is number of samples and m number of data(e.g sites)
    assumes data_filename format is

            sample_0, .., sample_n
    data_0
    .
    .
    data_m

    transpose before sending to function if you have different format

    missing_value_indicator - the missing value char (int, float or string) in your data 
    data_max_missing_values - the maximum data missing values allowed  (percantage - values between 0 and 1)
    samples_max_missing_values - the maximum sample missing values allowed  (percantage - values between 0 and 1)
    sep - the separator sign of the matrix in data_filename
    suffix - the suffix for the output filename

    returns array of non-missing data/samples
    """

    dim = 2

    #convert data_type from string to type ('float' --> float)
    data_type = DATA_TYPE
    original_data_type = DATA_TYPE

    
    #find the right missing value indicator type 
    replace = False
    float_ind = None
    int_ind = None
    try:
        # a string representing an int can be converted to both float and int
        # but a string representing a float cant be converted to int
        # that's why we check the float conversion first
        float_ind = float(missing_value_indicator) 
        int_ind = int(missing_value_indicator)
    except:
        pass
    # if both int_ind and float_int are not None - the string represent an int.
    # (if it represent an float -  it cant be converted to int) so its enough to check the int alone

    # if indicator type is not str - we'll use the original datatyp
    # if it is str - we'll read the data as strings, find the indicator and convert it back to datatype
    if int_ind is not None: 
        missing_value_indicator = int_ind
    elif float_ind is not None:
        missing_value_indicator = float_ind
    else:
        data_type = STR_DATA_TYPE
        replace = True

    try:
        all_data, col_names, row_names = common.load_data_file(data_filename, 2, header=header, na_values= missing_value_indicator)

        if all_data.ndim != dim:
            raw_input("Error: got data from dimensions %d while excepted to %d. Please check all paramenters are OK (data type and separator)." %(all_data.ndim, dim))
            return None
        output_filename = data_filename + suffix

        # if replace:
        #     data_type = original_data_type

        output_data, col_names, row_names = _replace_missing_values_in_matrix(all_data, missing_value_indicator, data_max_missing_values, samples_max_missing_values, replace, col_names, row_names)
        print "Output is saved to " + output_filename

        data_to_save = output_data
        if row_names is not None:
            data_to_save = column_stack((row_names, data_to_save))
        if col_names is not None:
            if row_names is not None:
                new_header = ["ID"] + list(col_names)
            else:
                new_header = col_names
            data_to_save = vstack((new_header, data_to_save))
            
        savetxt(output_filename, data_to_save, fmt='%s')
        return output_data

    except Exception as e:
        raw_input("Error loading data file. please check that data_type, separator and missing value indicator Ok\n%s"%e)
    

def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    

    parser.add_argument("--datafile", required=True, type=str, help="the matrix data filename")
    parser.add_argument("--chr", type=str, required=True, help="the missing value char (or any sign) in your data" )
    parser.add_argument("--maxs",required=True, type = float, help="the maximum missing values allowed per site (percentage - values between 0 and 1). If a site has more than this amount of missing values (samples) it will be deleted otherwise the missing values will be replaced by the mean of the site")
    parser.add_argument("--maxi", required=True,type = float, help="the maximum missing values allowed per sample (percentage - values between 0 and 1). If a sample has more than this amount of missing values (sites) it will be deleted.")
    parser.add_argument("--sep", type=str, default="\t", help="the separator sign of the matrix in data_filename. Default is tab")
    parser.add_argument("--suffix", type=str, default=".no_missing_values", help="the suffix for the output filename")
    # parser.add_argument("--header", action='store_true', help="supply this flag if your datafile contains a header")
    #parser.add_argument("--dim", type=int, default=2, help="the dimensions of the matrix in the datafile. Default is 2")
    

    return parser.parse_args(args=argv)

if __name__=="__main__":
    args = parse_args()
    replace_missing(args.datafile, args.chr, args.maxs, args.maxi, args.sep , args.suffix, None)
