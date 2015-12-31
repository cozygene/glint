import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname( os.path.abspath(__file__) ),"..",".."))
from modules import methylation_data
meth = methylation_data.MethylationData(datafile = args.datafile, excludefile = "datafile/excludefile.txt")
assertequal(meth, readmatrix( "datafile/excludefile.txt-res"))


def readmatrix(filepath) :
    from numpy import loadtxt
    data = loadtxt(filepath, dtype = str)
    samples_ids = data[0,:][1:]         # extract samples ID
    cpgnames = data[:,0][1:]           # extract methylation sites names
    
    return  data[1:,1:].astype(float) 