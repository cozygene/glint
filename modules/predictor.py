from numpy import matrix, array, column_stack, loadtxt,bincount, mean
from module import Module
class MethSite(object)
class Predictor(Module):
    def __init__(self):

    def handle_missing_values(data):
        masked_data = masked_array(data, isnan(data)) 
        mean_per_site = average(masked_data, axis=1)  
        # TODO is masked_data.mask equal to nan_indices? if so, we don't need to run this "where" line and just use masked_data.mask instead of nan_indices
        nan_indices = where(masked_data.mask) #TODO add [0] ?                 # find nan values indices
        self.data[nan_indices] = mean_per_site[nan_indices[0]]    # replace nan values by the mean of each site


    def bla():
            
        snps_info = load # [(site_id, allele, coeff),...]
        sites_info = load #[(site_score, number_of_predictors, site_snps_list), ...]

        # # follow left number of predictors to find # TODO check if could be 0
        # site_left_num_of_predictors = dict()
        # [site_left_num_of_predictors.setdefault(site_id, sites_to_snps[site_id][1]) for site_id in sites_to_snps.keys()]

        # # set initial prediction to 0
        # site_prediction = dict()
        # [site_prediction.setdefault(site_id,0) for site_id in sites_to_snps.keys()]

        # full_site_prediction = dict()

        # find only interesting methylation sites: sites who passed minimum score and those that we have all snps to predict
        plink_snps_ids = loadtxt(".map", dtype = str, usecols=(1))
        predicting_snps = set()
        interesting_sites = []
        for site_id in sites_info.keys():
            site_score = sites_info[site_id][0]
            site_snps = sites_info[site_id][2]
            # remove sites which score lower than minimum
            if site_score < min_score:
                sites_info.pop(site_id)
            # remove sites that we don't have the information about their predicting snps
            elif len(set(site_snps).difference(set(plink_snps_ids))) != 0:
                sites_info.pop(site_id)
            # if site is predictable - add it's snp to the interesting snp list
            else:
                interesting_sites.append(site_id)
                interesting_snps.update(site_snps)

        # find interesting snp columns indices in plink file
        interesting_snps = list(interesting_snps)
        snp_indices = []
        for i, snp in enumerate(interesting_snps):
            snp_index_in_plink_map = plink_snps_ids.index(snp)
            snp_indices.append(6 + snp_index_in_plink_map*2)
            snp_indices.append(6 + snp_index_in_plink_map*2 + 1)


        plink_data = loadtxt(".map", dtype = str, usecols=tuple(snp_indices.insert(0,1)))
        samples_list = plink_data[:,0]# TODO check
        plink_data = plink_data[:,1:]
        
        number_of_samples = len(samples_list) # TODO check
        

        meth_pred_columns_per_site = dict()
        [meth_pred_columns_per_site.setdefault(site, []) for site in interesting_sites]

        # calculate prediction
        for i in range(len(interesting_snps)):
            snp_info = snps_info[interesting_snps[i]]
            for site_id, allele, coeff in snp_info:
                samples_alleles = plink_data[:, [i, i+1]]
                snp_occurrences = bincount(where(samples_alleles==allele)[0], minlength = number_of_samples)
                
                # handle missing values - for samples with 2 nans : replace with std, sample with one nan (one nan allele) , keep 1
                nan_samples_indices = where(bincount(where(samples_alleles == 'NA')[0], minlength = number_of_samples) == 2) #samples who has 2 NAs
                not_nan_sample_indices = delete(range(number_of_samples), nan_samples_indices)
                # find avg accross not nans
                snp_occurrences_avg = mean(snp_occurrences[not_nan_sample_indices])
                # replace nans with avg
                snp_occurrences[nan_samples_indices] = snp_occurrences_avg
                # multiple by coeff
                predictor = coeff * snp_occurrences
                meth_pred_columns_per_site[site_id].append(predictor)

        
            
        predicted_meth_data = matrix([column_stack(tuple(meth_pred_columns_per_site[site_id])) for site_id in interesting_sites])
        meth_data_obj = MethylationData(predicted_meth_data, samples_list, interesting_sites, predicted = True)
        meth_data_obj.save("predicted_meth_data.glint")
        # TODO check that meth_data_obj.add_covar works. just run glint.py with the glintfile and add --covar

