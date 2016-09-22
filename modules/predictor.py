from numpy import vstack, loadtxt, bincount, mean, where, delete, zeros, concatenate, dtype, fromstring, seterr
import time
import logging
import pickle
import re
from module import Module
from methylation_data import MethylationData
from utils import common

class Predictor(Module):
    """
    predict methylation level of methylation sites by snps.
    the snps data given by plink files

    ignore sample (dont predict its methylation sites) when it has more than x snps that are missing (missing values)
    ignore snp (dont predict with it) when it has more than x samples that are missing (missing values)
    x has default value but can be changed
    
    ignore sites with low score (user can select the score) (score list is in the file "sites_scores_list")

    predict methylation level of a site if we have information about at least one snp that predicts it

    prediction model:
    If the model for some methylation site m is m=c1*s1+c2*s2 (model is given in the parsed data returned from convert_model_file_to_bin.py script
    a detailed explaintion is provided in PredictorParser)
    where c1,c2 are the coefficients in the model file and s1,s2 are the values of the SNPs s1,s2,
    then the predicted value of m for individual i will be m_i=c1*s1_i+c2*s2_i, where s1_i,s2_i are the values of SNPs s1,s2 in individual i
    """
    
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

        self.sites_name_per_id = loadtxt(sites_ids_list_file, dtype = str)

        # load snps coeffs for each site
        with open(sites_snps_coeff_list_file, 'r') as f:
            self.sites_snps_coeff = [[float(coeff) for coeff in line.split('\t')[:-1]] for line in f.read().splitlines()]
        f.close()

        self.sites_scores = loadtxt(sites_scores_list_file)

    def predict(self, min_score, plink_snp_file, plink_geno_file, plink_ind_file, min_missing_values):
        """
        predict with following rules:
        replace missing values with mean (unless there are more htan min_missing_values missing values)
        remove samples (dont predict them) which have more than min_missing_values missing snps (out of all its snps)
        remove snps (sont predict with them) which have more than min_missing_values missing samples (out of all its samples)
        remore site with score lower than min_score
        dont predict sites that we dont have any snp (if we have at least one - predict it)
        """
        samples = loadtxt(plink_ind_file, dtype=str, usecols=(0,))
        number_of_samples = samples.shape[0]
        if type(plink_ind_file) == file:
            plink_ind_file.close()

        plink_snps_data = loadtxt(plink_snp_file, dtype = str)
        if type(plink_snp_file) == file:
            plink_snp_file.close()
        # use only snps that their allele are not the pairs CG or AT - since in those cases we cannot know which strand was tested 
        relevant_snps_indices = self.get_relevant_plink_snp_list(plink_snps_data)

        logging.info("get number of snp occurences...")
        # extract find occurences per sample for each snp from .geno file
        snp_occurrences, relevant_snps_indices, missing_sampels_indices, non_missing_sampels_indices = self.get_snps_occurences(plink_geno_file, relevant_snps_indices, number_of_samples, min_missing_values) 
        
        # indices = [i  for i,name in enumerate(plink_snps_data[relevant_snps_indices,0]) if name in self.snps_id_per_name]
        relevant_snps_names = []
        relevant_snp_occurrences = []
        for i,name in enumerate(plink_snps_data[relevant_snps_indices,0]):
            if name in self.snps_id_per_name:
                relevant_snps_names.append(name)
                relevant_snp_occurrences.append(snp_occurrences[i])

        # remove snps that we dont have information on
        self.predicted_samples = samples[non_missing_sampels_indices]
        number_of_samples = self.predicted_samples.size

        if (number_of_samples == 0):
            common.terminate("All samples removed. There is nothing to predict. quiting...")

        # find sites with score bigger than min_score
        seterr(invalid='ignore') # to ignore the following line warning (These warnings are an intentional aspect of numpy)
        relevant_sites_indices = where(self.sites_scores > min_score)[0]
        logging.info("remove %s sites with score lower than %s..." % (len(relevant_sites_indices), min_score))

        # this code removes the sites in the list bad_sites_list (list of cpgs) from our data. it wast tested.
        # logging.info("remove bad sites..." )
        # bad_sites_indices = where(in1d(self.sites_name_per_id, bad_sites_list))[0]
        # relevant_sites_indices = delete(relevant_sites_indices, bad_sites_indices)

        # calc prediction
        logging.info("predict methylation level...")
        site_prediction, predicted_sites_ids =  self.predict_sites(number_of_samples, relevant_snps_names, relevant_snp_occurrences, relevant_sites_indices)
        if site_prediction == []: # no sites predicted
            common.terminate("All sites removed. There is nothing to predict.")
        
        self.predicted_sites_names = self.sites_name_per_id[predicted_sites_ids]
        self.site_prediction = site_prediction
        
    def get_relevant_plink_snp_list(self, plink_snps_data):
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

    def get_snps_occurences(self, plink_geno_file, snp_indices, number_of_samples, min_missing_values):
        """
        extracts number of occerences for each snp and each sample from .geno plink file, handles missing values
        the format of this file:
            there is a line for each SNP
            the line is a string above [0,1,2] which represent number of occurences of this snp for each sample 
            (number of occurences of snp i in sample j is the value in line i in the j'th index)

        if there are too many missing values in snp X: ignore this snp for all samples, and dont use it to predict site methylation level
        otherwise, predict it's missing values by average over all other samples number of occurences 

        if there are too  many missing values is sample X: don't predict this sample methylation levels at all

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
        missing_sampels_indices = where(na_count_per_sample > number_of_samples * min_missing_values)[0]
        non_missing_sampels_indices = delete(range(number_of_samples), missing_sampels_indices)

        # remove missing samples from snp occurences
        relevant_snp_occurrences = vstack(tuple(relevant_snp_occurrences))
        relevant_snp_occurrences = relevant_snp_occurrences[:, non_missing_sampels_indices]

        #missing_sampels_indices = where(na_count_per_sample >= number_of_samples * 0)[0]
        logging.info("removing %d samples with more than %f missing values..." %(len(missing_sampels_indices), min_missing_values))
        
        return relevant_snp_occurrences, relevant_snps_indices, missing_sampels_indices, non_missing_sampels_indices


    def predict_site(self, number_of_samples, site_snps_ids, snps_coeffs, relevant_snps_ids, relevant_snp_occurrences):
        """
        site_snps_ids - the ids of all snps that are the site predictors
        snps_coeffs - an array such that snps_coeffs[i] is the coefficient of the snp which id is site_snps_ids[i]
        relevant_snps_ids - list of the relevant snps for prediction (we don't predict using 
                            every snp, for wxample we ignore snps which alleles are A and T)
        relevant_snp_occurrences - an array of arrays , size aXb where a is number of relevant snps and b is number of samples.
                                    relevant_snp_occurrences[i][j] is number of alleles of 
                                    snp id ID such that relevant_snps_indices_per_id[ID] = i  in sample j.

        predict each site with the following model:

        Assume we have a methylation site m with two predictors s1, s2 and assume that their reference alleles are G, C and coefficients c1, c2,
        respectively. Given the genotypes of an individual i in SNPs s1,s2  (denote s1_i,s2_i), we predict m_i, the methylation level of 
        individual i in meth site m, as follows:
        m_i_predicted = c1*(number of 'G' alleles in s1_i) + c2*(number of 'C' alleles in s2_i)
        """
        snp_index_per_id = dict()
        for index, snp_id in enumerate(relevant_snps_ids):
            snp_index_per_id[snp_id] = index

        site_prediction = zeros(number_of_samples)

        for i,snp_id in enumerate(site_snps_ids):
            if snp_id in relevant_snps_ids:
                site_prediction += snps_coeffs[i] * relevant_snp_occurrences[snp_index_per_id[snp_id]]

        return site_prediction


    def predict_sites(self, number_of_samples, relevant_snps_names, relevant_snp_occurrences, relevant_sites_indices):

        """
        predict each site that we have  at least informations about one of the sites predicting snps
        if we have no snp to predict with - dont predict the site

        number_of_samples - (int) number of samples (n)
        relevant_snps_names - an array of names of the relevant snps for sites prediction (we don't predict using 
                              every snp, for example we ignore snps which alleles are A and T).
                              array of size s - number of relevant snps
                              relevant_snps_names[i] is the name of the ith snp ("rsXXX..")
        relevant_snp_occurrences - array of size sXn where s is the number of snps that are relevant for this site prediction
                                    and  n is number of samples.
                                    relevant_snp_occurrences[i][j] is the number of occurences of the i-th snp in sample j.
        relevant_sites_indices - a list (array) of the indices of the sites that we will predict now. array of size m - number of sites

        return:
            matrix of size mXn (m- number of predicted sites, n- number of samples) of the prediction
            predicted_sites_ids - list of the predicted sited ids

        """
        
        relevant_snps_ids = [self.snps_id_per_name[name] for name in relevant_snps_names]
        relevant_snps_ids_set = set(relevant_snps_ids)

        # iterate over each site and predict its value
        next_index = 0
        sites_predictions = []
        predicted_sites_ids = []

        with open(self.site_snps_list_file, 'r') as f:
            # iterate over only relevant sites
            for i, site_snps in enumerate(f):

                if (next_index < len(relevant_sites_indices)) and (i == relevant_sites_indices[next_index]):
                    site_snps_ids = [int(sid) for sid in site_snps.split('\t')[:-1]]
                    
                    # if we have information about at least one of the sites predicting snps
                    if len(set(site_snps_ids).difference(relevant_snps_ids_set)) < len(site_snps_ids):
                        site_prediction = self.predict_site(number_of_samples, site_snps_ids, self.sites_snps_coeff[i], relevant_snps_ids, relevant_snp_occurrences)
                        sites_predictions.append(site_prediction)
                        predicted_sites_ids.append(i)

                    next_index += 1

        f.close()
        if len(sites_predictions) == 0:
            return [], []
        return vstack(tuple(sites_predictions)), predicted_sites_ids


    def meth_data(self):
        return MethylationData(self.site_prediction, self.predicted_samples, self.predicted_sites_names)
        



