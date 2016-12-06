import os
from numpy import where, vstack, savetxt
from pandas import Index, unique, read_csv, DataFrame

"""
this script creates the HumanMethylationSites file (you can choose output file name with outputfile parameter)
it also saves a list of cpg sites found in X and Y chromosomes in a seperate file

it extracts the fields from the chip data
it first reads the data from 850K than only adds information about new cpg sites from the rest of the files
(it takes all the sites from file files[0] than adds the different sites between files[1] and files[0],
    than adds only the different sites between files[2] and files[1]....)

The first file in the files list is assumed to be the most updated one (the one to first take the information from)

to add new chip, add another filepath (like filepath_850) and add it to the list of files ("files" parameter) in the right location 
where the first file (location 0) is considered the most updated data file
"""

t=['IlmnID','CHR','MAPINFO','UCSC_RefGene_Name','Relation_to_UCSC_CpG_Island'] # the interesting fields
repeat_field = 'UCSC_RefGene_Name'

#paths to the chip dataset
filepath_850 = '/Users/yedidimr/Downloads/MethylationEPIC_v-1-0_B2.csv'
filepath_450 = '/Users/yedidimr/Downloads/HumanMethylation450_15017482_v1-2.csv'
filepath_270 = '/Users/yedidimr/Downloads/humanmethylation27_270596_v1-2.bpm'
files = [filepath_850, filepath_450]

# path to the output file
outputfile = "HumanMethylationSites_new"


def get_data(filepath):
    """
    reads a chip 
    """
    f=open(filepath,'r')
    i=0
    while(i<8):
         i+=1
         dd = f.readline().split(",")
    f.close()

    # first 7 rows are headers, skip them
    # read only the columns of interest
    data = DataFrame.as_matrix(read_csv(filepath, dtype = str, sep=',', usecols=[dd.index(r) for r in t], skiprows = 7))    
    data = data.astype(str)
    return data

def new_cpgs(old_data, new_data):
    """
    cpg are assumed to be in the first column
    return the indexes of the cpg ids in the new_data that are not found in the old_data
    """
    old_cpgid =  old_data[:,0]
    new_cpgid =  new_data[:,0]

    # find indexes cpg sites in the new data which do not exist in the old data
    indexes_list = where(Index(unique(old_cpgid)).get_indexer(new_cpgid) <0)[0]
    return indexes_list


"""
read all the files
"""
data = None
for i,datafile in enumerate(files):
    print "loading %s "% datafile

    if data is None:
        data = get_data(datafile)
    else:
        new_data = get_data(datafile)
        
        # find indexes cpg sites in 450K which do not exist in 850K
        new_cpg_indexes =  new_cpgs(data, new_data)
        if len(new_cpg_indexes) > 0:
            print "found %s cpgs in %s which do not exist in the more updated datafiles (those who were previous loaded)" % (len(new_cpg_indexes) , os.path.basename(datafile))
            data = vstack((data, new_data[new_cpg_indexes,:]))
            

"""
remove repeats in the UCSC field:
for example in cpg id cg00214611 the value is "TMSB4Y;TMSB4Y" so we'll change it to "TMSB4Y"
"""
repeat_i = t.index(repeat_field)
for i in range(data.shape[0]):
    data[i, repeat_i] = ";".join(list(set(data[i,repeat_i].split(";"))))

print "saving output to %s" % outputfile

savetxt(outputfile, data, delimiter=",",header = ",".join(t), fmt="%s")

# save X Y chromosomes
y_chr_ind = list(where(data[:,1]=='Y')[0])
x_chr_ind = list(where(data[:,1]=='X')[0])
xy_ind = y_chr_ind
xy_ind.extend(x_chr_ind)
xy_chr_output_file = outputfile +"_XY_chr_only"
print "Saving cpg sites for X Y chromosomes only to %s" % xy_chr_output_file
savetxt(xy_chr_output_file, data[xy_ind, 0].reshape(-1,1), fmt="%s")


