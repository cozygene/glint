from numpy import vstack, loadtxt, bincount, mean, where, delete, zeros, concatenate, dtype, fromstring
import time
import logging
import pickle
import re
from module import Module
from methylation_data import MethylationData


class Predictor(Module):
    NA_VALUE = 9
    SITES_SCORES_FILE = 'sites_scores_list'
    SITES_SNPS_FILE = 'site_snps_list'
    SITES_IDS_FILE ='sites_ids_list'
    SNPS_IDS_FILE = 'snps_ids_list'
    SITES_SNPS_COEFF_FILE = 'sites_snps_coeff_list'
    
    def __init__(self, sites_scores_list_file, site_snps_list_file,  sites_ids_list_file, snps_ids_list_file, sites_snps_coeff_list_file):

        logging.info("loading model files...")
        self.site_snps_list_file = site_snps_list_file # will load it line by line later 

        with open(snps_ids_list_file, 'r') as f:
            self.snps_id_per_name = dict()
            for snp_id, snp_name in enumerate(f.read().splitlines()):
                self.snps_id_per_name[snp_name] = snp_id
        f.close()

        self.sites_id_per_name = loadtxt(sites_ids_list_file, dtype = str)

        # load snps coeffs for each site
        with open(sites_snps_coeff_list_file, 'r') as f:
            self.sites_snps_coeff = [[float(coeff) for coeff in line.split('\t')[:-1]] for line in f.read().splitlines()]
        f.close()

        self.sites_scores = loadtxt(sites_scores_list_file)

    def predict(self, min_score, plink_snp_file, plink_geno_file, plink_ind_file, min_missing_values):
        samples = loadtxt(plink_ind_file, dtype=str, usecols=(0,))
        number_of_samples = samples.shape[0]
        plink_ind_file.close()

        plink_snps_data = loadtxt(plink_snp_file, dtype = str)
        plink_snp_file.close()
        # snps that dont containg CG or AT
        relevant_snps_indices = self.get_plink_snp_list(plink_snps_data)

        logging.info("get number of snp occurences...")
        # extract find occurences per sample for each snp from .geno file
        relevant_snp_occurrences, relevant_snps_indices, missing_sampels_indices = self.get_relevant_snps_occurences(plink_geno_file, relevant_snps_indices, number_of_samples, min_missing_values) 
        relevant_snps_names = plink_snps_data[relevant_snps_indices, 0]
        non_missing_sampels_indices = delete(range(number_of_samples), missing_sampels_indices)

        # number_of_samples = len(non_missing_sampels_indices)
        self.predicted_samples = samples[non_missing_sampels_indices]
        # find sites with score bigger than min_score
        logging.info("remove sites with score lower than %s..." % min_score)

        relevant_sites_indices = where(self.sites_scores > min_score)[0]

        # calc prediction
        logging.info("predict methylation level...")
        site_prediction, predicted_sites_ids =  self.predict_sites(number_of_samples, relevant_snps_names, relevant_snp_occurrences, relevant_sites_indices, relevant_snps_indices)
        self.predicted_sites_names = self.sites_id_per_name[predicted_sites_ids]
        logging.info("remove samples with more than %s missing values..." % min_missing_values)
        
        # remove missing samples from prediction
        self.site_prediction = site_prediction[:, non_missing_sampels_indices]
        

    def get_plink_snp_list(self, plink_snps_data):
        """
        gets a .snp file and deletes all snps which references (two last columns in the file) are ('C' and 'G') or ('A' and 'T')
        the snps that will be left are those which references are ('C' and 'A') or ('C' and 'T') or ('G' and 'A') or ('G' and 'T')
        we want to delete those snps since we cann't know which strand was tested

        plink_snps_data - a .snp eigenstrat plink data

        returns the indicis of the left snps
        """
        
        num_of_snps = plink_snps_data.shape[0]

        c = bincount(where(plink_snps_data == 'C')[0], minlength=num_of_snps)
        a = bincount(where(plink_snps_data == 'A')[0], minlength=num_of_snps)
        t = bincount(where(plink_snps_data == 'T')[0], minlength=num_of_snps)
        g = bincount(where(plink_snps_data == 'G')[0], minlength=num_of_snps)

        cg_indices = where(c&g == 1)[0]
        at_indices = where(a&t == 1)[0]

        indices_to_delete = concatenate((at_indices, cg_indices), axis=0)
        relevant_indices = delete(range(num_of_snps), indices_to_delete)
        return relevant_indices


    def convert_012_string_to_ndarray(self, string):
        return (fromstring(string, dtype=dtype('S1'))).astype(int)

    def get_relevant_snps_occurences(self, plink_geno_file, snp_indices, number_of_samples, min_missing_values):
        """
        extracts number of occerences for each snp and each sample from .geno plink file, handles missing values
        the format of this file:
            there is a line for each SNP
            the line is a string above [0,1,2] which represent number of occurences of this snp for each sample 
            (number of occurences of snp i in sample j is the value in line i in the j'th index)

        if there are too many missing values in snp X: ignore this snp for all samples, and dont use it to predict site methylation level
        otherwise, predict it's missing values by average over all other samples number of occurences 

        if there are too  many missing values is sample X: don't predict this sample methylation levels at all

        TODO -if there are too many snp that we ignored in a site?

        input:
        plink_geno_file - path to the .geno plink file
        snp_indices - list of snp indices that we should extract from the file, each index represent a line in the .geno file (the i'th value in snp_indices is the snp_indices[i]'th line in .geno file)

        output:
        relevant_snp_occurrences - a list of length of number of relevant snps that we used (snp_indices without too missing snps that we found in this function)
                                   each cell represent a SNP and contains a array of number of occurences of this snp in each sample (extracted from .geno file)
        relevant_snps_indices - list of the snp indices that we extracted 
        missing_sampels_indices - list of too missing samples indices

        """
        relevant_snp_occurrences = [] # will hold ndarray of occurences per sample for each snp
        na_count_per_sample  = zeros(number_of_samples)

        if type(plink_geno_file) != file:
            samples_snps_f = open(plink_geno_file, 'r')
        else:
            samples_snps_f = plink_geno_file

        snps_missing_values_counter = 0
        samples_missing_values_counter = 0

        # iterate over each snp
        next_index = 0
        for i, samples_snp in enumerate(samples_snps_f):
            if (next_index < len(snp_indices)) and (i == snp_indices[next_index]):

                snp_occurrences = self.convert_012_string_to_ndarray(samples_snp[:-1])
                na_indices = snp_occurrences[where(snp_occurrences == self.NA_VALUE)[0]] #samples indices where this snp is missing
                na_count_per_sample[na_indices] += 1 

                na_count = len(na_indices)
                na_percentage = float(na_count) / number_of_samples
                if na_percentage > min_missing_values:
                    # too many missing values in this snp: ignore this snp for all samples, and dont use it to predict site methylation level
                    snp_indices[next_index] = -1 # mark as not relevant
                    snps_missing_values_counter += 1
                else:
                    if na_count != 0 :
                        # "predict" snp occurences  - relate the mean of non-missing samples in this snp
                        non_na_indices = delete(range(len(snp_occurrences)), na_indices)
                        snp_occurrences[na_indices] = mean(snp_occurrences[non_na_indices])

                    relevant_snp_occurrences.append(snp_occurrences)

                next_index += 1

        samples_snps_f.close()
        
        logging.info("removing %d snps with more than %f missing values..." %(snps_missing_values_counter, min_missing_values))
        relevant_snps_indices = snp_indices[where(snp_indices != -1)[0]]
        
        # find samples with too many missing values
        missing_sampels_indices = where(na_count_per_sample > number_of_samples * min_missing_values)
        logging.info("removing %d samples with more than %f missing values..." %(len(missing_sampels_indices), min_missing_values))
        
        return relevant_snp_occurrences, relevant_snps_indices, missing_sampels_indices

    def predict_site(self, number_of_samples, site_snps_ids, snps_coeffs, relevant_snps_indices_per_id, relevant_snp_occurrences):
        site_prediction = zeros(number_of_samples)

        for j, snp_id in enumerate(site_snps_ids):
            try:
                snp_index = relevant_snps_indices_per_id[snp_id]
            except:
                continue
            site_prediction += snps_coeffs[j] * relevant_snp_occurrences[snp_index]
        return site_prediction


    def predict_sites(self, number_of_samples, relevant_snps_names, relevant_snp_occurrences, relevant_sites_indices, relevant_snps_indices):
        # look at sites that we have all snps to predict
        relevant_snps_indices_per_id = dict()

        relevant_snps_number = len(relevant_snps_indices)
        for i in range(relevant_snps_number):
            try:
                snp_id = self.snps_id_per_name[relevant_snps_names[i]] 
            except:
                continue # plink snp not in our model
            relevant_snps_indices_per_id[snp_id] = i

        relevant_snps_ids_set = set(relevant_snps_indices_per_id.keys())

        next_index = 0
        sites_predictions = []
        predicted_sites_ids = []
        with open(self.site_snps_list_file, 'r') as f:
            for i, site_snps in enumerate(f):
                                
                if (next_index < len(relevant_sites_indices)) and (i == relevant_sites_indices[next_index]):
                    site_snps_ids = [int(sid) for sid in site_snps.split('\t')[:-1]]
                    
                    if len(set(site_snps_ids).difference(relevant_snps_ids_set)) < len(site_snps_ids):

                        site_prediction = self.predict_site(number_of_samples, site_snps_ids, self.sites_snps_coeff[i], relevant_snps_indices_per_id, relevant_snp_occurrences)
                        sites_predictions.append(site_prediction)
                        predicted_sites_ids.append(i)

                    next_index += 1
        f.close()
        return vstack(tuple(sites_predictions)), predicted_sites_ids


    def meth_data(self):
        return MethylationData(self.site_prediction, self.predicted_samples, self.predicted_sites_names)
        



