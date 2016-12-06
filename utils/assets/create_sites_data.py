from numpy import where, vstack, savetxt
from pandas import Index, unique
from pandas import read_csv, DataFrame

t=['IlmnID','CHR','MAPINFO','UCSC_RefGene_Name','Relation_to_UCSC_CpG_Island'] # the interesting fields
repeat_field = 'UCSC_RefGene_Name'

#paths to the chip dataset
filepath_850 = '/Users/yedidimr/Downloads/MethylationEPIC_v-1-0_B2.csv'
filepath_450 = '/Users/yedidimr/Downloads/HumanMethylation450_15017482_v1-2.csv'
filepath_270 = '/Users/yedidimr/Downloads/humanmethylation27_270596_v1-2.bpm'

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


print "loading %s "% filepath_850
data850 = get_data(filepath_850)
cpgid850 =  data850[:,0]

print "loading %s "% filepath_450
data450 = get_data(filepath_450)
cpgid450 =  data450[:,0]

# find indexes cpg sites in 450K which do not exist in 850K
indexes_list = where(Index(unique(cpgid850)).get_indexer(cpgid450) <0)[0]
print "found %s cpgs in 450K which do not exist in 850K" % len(indexes_list) 
if len(indexes_list) != 0: 
    data = vstack((data850, data450[indexes_list,:]))
else:
    data = data850

"""
remove repeats in the UCSC field:
for example in cpg id cg00214611 the value is "TMSB4Y;TMSB4Y" so we'll change it to "TMSB4Y"
"""
repeat_i = t.index(repeat_field)
for i in range(data.shape[0]):
    data[i, repeat_i] = ";".join(list(set(data[i,repeat_i].split(";"))))

print "saving output to %s" % outputfile

savetxt(outputfile, data, delimiter=",",header = ",".join(t), fmt="%s")