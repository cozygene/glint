from numpy import array, delete, where, in1d

"""
This script is not for glint users but for the SW developers.
It does the preprocessing for imputation of methylation level by SNPs.

This script parsers the "model"* file, i.e. the file with the coefficients for making methylation imputations from SNPs.
The file is "heavy" and it takes time to read and parse it. So this preprocessing-script does the "hard work",
extracts and saves the relevant information in a way that it will be efficient to get it on "real time".

The output files of the deleted scripts are
    - "sites_ids_list" a list of cpgs (that can be explained by SNPs), which describes the cpgs IDs: the id of the cpg name at line number i is i.
    - "snps_ids_list"  a list of SNPs descibing their ids:  the id of the SNP at line number i is i.
    - "site_snps_list" which describe for each site the snps that "explain" it. the numbers in the i'th line are the ids of the snps descibe the site which id is i.
    - "sites_scores_list" the score at the i'th line is the imputations correlation of i'th site (site which id is i)
    - "sites_snps_coeff_list" the numbers at the i'th line are the coefficients of the snps which "explain" the i'th site. the j'th number is the coefficient of the j'th snp at the i'th i in the file "site_snps_list"

Note: I decided to save lists and not pickles after comparing the runtime of loading both of them. 
"loadtxt" is used for loading lists, and it makes the loading much faster than pickle.

This script should be executed once, in cases when:
- if you have different model file than the one we used for the last version.
- if you need to extract more data than we already extracted.

* in this case the "model" file name was KORA_model_multiple_snps_W_50_M_10.
 The structure of the file, culumn by column: 
  cpg id, chromosome num, score (imputation correlation), lasso's lambda, number of imputing SNPs (an integer in the range 1-10),
  the position in basepairs of the imputing SNPs (the average across the imputing SNPs).
  The following columns are triplets, each describing one of the the imputing SNPs of the current cpg. 
  For each imputing SNP we have 3 values: "rs" id (the identifier of the SNP - note that it doesn't always start with "rs"),
  reference allele, coefficient.
"""
# the path to the model file which will be preproccesed 
KORA_MODEL_FILE_PATH =  "enter the path here" #KORA_model_multiple_snps_W_50_M_10.txt" 

 # path to the list artifacts sites
BAD_PROBES_FILE_PATH = "parsers/assets/polymorphic_cpgs.txt"
SITE_ID_INDEX = 0
IMPUTATION_CORRELATION_SCORE_INDEX = 2
NUMBER_OF_IMPUTATING_SNPS_INDEX = 4
FIRST_SNP_INDEX = 6

data = file(KORA_MODEL_FILE_PATH, 'r').read()
artifacts = array(file(BAD_PROBES_FILE_PATH, 'r').read().splitlines())

sites = data.splitlines()
sites_names = array([site.split('\t')[0] for site in sites])

# remove artifacts
sites_indices_to_remove = where(in1d(sites_names , artifacts))[0]

sites = array(sites)
keep_sites = delete(sites, sites_indices_to_remove)

sites_ids = []
snps_ids = []
sites_scores = []
current_site_snps_list = []
sites_snps_list = []
current_site_coeffs_list = []
sites_coeffs_list = []
current_snp_index = 0
print "Parsing data... this gonna take some time..."

for site in keep_sites:
    site_info = site.split('\t')
    site_id = site_info[SITE_ID_INDEX]
    site_score = site_info[IMPUTATION_CORRELATION_SCORE_INDEX]
    number_of_predictors = int(site_info[NUMBER_OF_IMPUTATING_SNPS_INDEX])
    current_site_snps_list = []
    current_site_coeffs_list = []
    for i in range(number_of_predictors):
        snp_id = site_info[FIRST_SNP_INDEX + i*3]
        coeff = float(site_info[FIRST_SNP_INDEX + i*3 + 2])
        
        if snp_id not in snps_ids:
            snps_ids.append(snp_id)
            snp_index = current_snp_index
            current_snp_index += 1
        else:
            snp_index = snps_ids.index(snp_id)

        current_site_snps_list.append(str(snp_index))
        current_site_coeffs_list.append(str(coeff))

    sites_ids.append(site_id)
    sites_scores.append(str(site_score))
    sites_snps_list.append("\t".join(current_site_snps_list))
    sites_coeffs_list.append("\t".join(current_site_coeffs_list))


with open("sites_snps_coeff_list", 'wb') as f1:
    f1.write("\n".join(sites_coeffs_list))
f1.close()

with open("site_snps_list", 'wb') as f1:
    f1.write("\n".join(sites_snps_list))
f1.close()


with open("sites_scores_list", 'wb') as f1:
    f1.write("\n".join(sites_scores))
f1.close()

with open("sites_ids_list", 'wb') as f1:
    f1.write("\n".join(sites_ids))
f1.close()

with open("snps_ids_list", 'wb') as f1:
    f1.write("\n".join(snps_ids))
f1.close()
