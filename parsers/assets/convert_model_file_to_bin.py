from numpy import array, delete
from pickle import dump



SITE_ID_INDEX = 0
PREDICTION_CORRELATION_SCORE_INDEX = 2
NUMBER_OF_PREDICTIONG_SNPS_INDEX = 4
FIRST_SNP_INDEX = 6
data = file("parsers/assets/KORA_model_multiple_snps_W_50_M_10.txt", 'r').read()
artifacts = file("parsers/assets/artifacts_chen.2013.txt", 'r').read().splitlines()
sites = data.splitlines()
sites_names = [site.split('\t')[0] for site in sites]
sites_indices_to_remove = [i for i in range(len(sites_names)) if sites_names[i] in artifacts]
sites = array(sites)
keep_sites = delete(sites, sites_indices_to_remove)


snps_info = dict()
sites_info = dict()
for site in keep_sites:
    site_info = site.split('\t')
    site_id = site_info[SITE_ID_INDEX]
    site_score = site_info[PREDICTION_CORRELATION_SCORE_INDEX]
    number_of_predictors = int(site_info[NUMBER_OF_PREDICTIONG_SNPS_INDEX])
    site_snps_list = []
    for i in range(number_of_predictors):
        snp_id = site_info[FIRST_SNP_INDEX + i*3]
        allele = site_info[FIRST_SNP_INDEX + i*3 + 1]
        coeff = float(site_info[FIRST_SNP_INDEX + i*3 + 2])
        site_snps_list.append(snp_id)
        if snp_id not in snps_info:
            snps_info[snp_id] = [(site_id, allele, coeff)]
        else:
            snps_info[snp_id].append((site_id, allele, coeff))

    sites_info[site_id] = (site_score, number_of_predictors, site_snps_list) #TODO remove unneccecery site_snps_list

with open("model-snps_info.pickle", 'wb') as f1:
    dump(snps_info, f1)
f1.close()


with open("model-sites_info.pickle", 'wb') as f2:
    dump(sites_info, f2)
f2.close()