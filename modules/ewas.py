import os
from utils import regression, tools#, plot
from numpy import column_stack, ones, savetxt, array, insert, vstack, loadtxt, append, where
from module import Module
from utils import common, plot, sitesinfo
import logging

"""
copy meth_data in advance
"""
class LinearRegression(Module):
    def __init__(self, methylation_data):
        self.meth_data = methylation_data

    def run(self):
        logging.info('running linear regression test...');
        #running association tests
        res =  self._linear_regression_test()
        logging.info('EWAS is Done!')
        return res

        
    def _linear_regression_test(self, output_filename = None):
        """
        linear regression test
        returns output sorted by p-values 

        return values:
            sorted_cpgnames - a list of cpgnames 
            sorted_pvalues - a list of p values 
            sorted_fstats - list of t-statistic 
            sorted_intercept_beta - coefficients of the intercept
            sorted_covars_betas - coefficients of the covariates
            sorted_site_beta - coefficients of the site under test

        sorted_pvalues[i] is the pvalue of the site sorted_cpgnames[i] (sorted_fstats[i] is its statistic and so on..)
        """
        logging.info("running linear regression test...")
        output = []           
        for i, site in enumerate(self.meth_data.data):
            coefs, fstats, p_value = regression.LinearRegression.fit_model(self.meth_data.phenotype, site, covars = self.meth_data.covar)
            
            # Note: if you add more info to site info note to:
            #       -   keep p_value at index 1 since the data is sorted by index 1
            #       -   increase / decrease number of values in the line if  output.shape[1] = X
            site_info = [self.meth_data.cpgnames[i], p_value[-1], fstats[-1]]
            site_info.extend([coefs[i] for i in range(coefs.size)]) 
            output.append(site_info) 
            

        
        output.sort(key = lambda x: x[1]) # sort output by p-value (1 is p-value index)
        output = array(output)
        sorted_cpgnames = output[:,0]
        sorted_pvalues  = output[:,1].astype(float)
        sorted_fstats   = output[:,2].astype(float)
        sorted_intercept_beta = output[:,3].astype(float)
        sorted_site_beta      = output[:,-1].astype(float)
        if  output.shape[1] == 5: # there is no covariates coefficient
            sorted_covars_betas = None
        else:
            sorted_covars_betas   = output[:,4:-1].astype(float)
        return sorted_cpgnames , sorted_pvalues, sorted_fstats, sorted_intercept_beta, sorted_covars_betas, sorted_site_beta


class EWASResults(object):
    CPGNAMES_INDEX = 0
    CHR_INDEX = 1
    POSITION_INDEX = 2
    PVALUE_INDEX = 3
    QVALUE_INDEX = 4
    GENE_INDEX = -2
    CATEGORY_INDEX = -1 #island

    CPGNAMES_TITLE = "{test_name}:ID"
    CHR_TITLE = "chromosome"
    POSITION_TITLE = "MAPINFO" #position
    PVALUE_TITLE = "p-values"
    QVALUE_TITLE = "q-values"
    GENE_TITLE = "UCSC_RefGene_Name" #gene
    CATEGORY_TITLE = "Relation_to_UCSC_CpG_Island" #island
    STATS_TITLE = "statistic" # statistics
    SIGMA_G_TITLE = "sigma-g"
    SIGMA_E_TITLE = "sigma-e"
    INTERCEPT_TITLE = "intercept" # the intercept coefficients
    SINGLE_COVAR_TITLE = "V{covar_index}" # the title for covar No1 is "V1"
    BETA_TITLE = "beta" # the site under test coefficients

    DELIMITER = ','

    def __init__(self, test_name, cpgnames, pvalues, qvalues = None, statistic = None, intercept_coefs = None, covars_coefs = None, site_coefs = None, sigma_g = None, sigma_e = None, sites_info_obj = None):
        
        self.test_name = test_name
        self.cpgnames = cpgnames
        self.pvalues = pvalues
        self.stats = statistic
        self.intercept = intercept_coefs
        self.beta = site_coefs
        self.sigma_g = sigma_g
        self.sigma_e = sigma_e

        self.covariates_betas = covars_coefs
        if covars_coefs is not None:
            if covars_coefs.ndim == 1:
                self.covariates_betas = covars_coefs.reshape((-1,1)) # make sure dim is (x,1) and not (x,)
            self.num_of_covars = covars_coefs.shape[1]
        else:
            self.num_of_covars = 0

        if qvalues is None:
            self.qvalues = tools.FDR(pvalues)
        else:
            self.qvalues = qvalues
        
        if sites_info_obj is None:   
            self.sites_info = sitesinfo.SitesInfoGenerator(self.cpgnames)
        else:
            self.sites_info = sites_info_obj

    def get_covars_titles(self, num_of_covars):
        return [self.SINGLE_COVAR_TITLE.format(covar_index = i+1) for i in range(num_of_covars)]



class EWASResultsCreator(EWASResults):
    def __init__(self, test_name, cpgnames, pvalues, qvalues = None, statistic = None, intercept_coefs = None,\
                 covars_coefs = None, site_coefs = None, sigma_g = None, sigma_e = None, sites_info_obj = None):
        """
        all inputs must be numpy arrays!

        cpgnames - array (numpy array) of cpgnames
        pvalues - array (numpy array) of pvalues for each cpg (i.e pvalues[i] is the pvalue of site named cpgnames[i])
        
        optional arguments:
        qvalues -  the qvalues array
        statistic - array with the t-statistic for each site
        intercept_coefs - array with the coefficients of the intercepts
        covars_coefs - array or matrix with the coefficients for each covariates
        site_coefs - the coefficients of the site under test
        sigma_g - lmm value for each site
        sigma_e - lmm value for each site
        sites_info_obj - SiteInfo bject describing the cpgnames
        """
        logging.info("Generating %s results file..." % test_name)
        super(EWASResultsCreator, self).__init__(test_name, cpgnames, pvalues, qvalues, statistic, intercept_coefs, covars_coefs, site_coefs, sigma_g, sigma_e, sites_info_obj)
        self.title = self.generate_title()
        self.data = self.generate_data()
                

    def generate_title(self):
        title = [self.CPGNAMES_TITLE.format(test_name = self.test_name), self.CHR_TITLE, self.POSITION_TITLE, self.PVALUE_TITLE, self.QVALUE_TITLE]

        if self.intercept is not None:
            title.append(self.INTERCEPT_TITLE)

        if self.covariates_betas is not None:
            covars_titles = self.get_covars_titles(self.num_of_covars)
            title.extend(covars_titles)

        if self.beta is not None:
            title.append(self.BETA_TITLE)

        if self.stats is not None:
            title.append(self.STATS_TITLE)

        if self.sigma_e is not None:
            title.append(self.SIGMA_E_TITLE)

        if self.sigma_g is not None:
            title.append(self.SIGMA_G_TITLE)

        title.extend([self.GENE_TITLE, self.CATEGORY_TITLE])

        return array(title)


    def generate_data(self):
        data_columns = [self.cpgnames, self.sites_info.chromosomes, self.sites_info.positions, self.pvalues, self.qvalues]

        if self.intercept is not None:
            data_columns.append(self.intercept)

        if self.covariates_betas is not None:
            [data_columns.append(self.covariates_betas[:,i]) for i in range(self.num_of_covars)] 

        if self.beta is not None:
            data_columns.append(self.beta)

        if self.stats is not None:
            data_columns.append(self.stats)

        if self.sigma_e is not None:
            data_columns.append(self.sigma_e)

        if self.sigma_g is not None:
            data_columns.append(self.sigma_g)

        data_columns.extend([self.sites_info.genes, self.sites_info.categories])
        return array(data_columns).T
        
    def save(self, output_filename):
        output = vstack((self.title, self.data))
        logging.info("%s results are saved to file %s" %(self.test_name, output_filename))
        savetxt(output_filename, output, fmt = "%s", delimiter=self.DELIMITER)

        
class EWASResultsParser(EWASResults):
    """docstring for EWASResultsParser"""
    def __init__(self, results_filemame):
        logging.info("Reading resultds from file %s" % results_filemame)
        data = self.readfile(results_filemame)
        test_name, cpgnames, pvalues, qvalues, stats, intercept, covariates_betas, beta, sigma_g , sigma_e, sitesinfo = self.parsedata(data)
        super(EWASResultsParser, self).__init__(test_name, cpgnames, pvalues, qvalues, stats, intercept, covariates_betas, beta, sigma_g , sigma_e, sitesinfo)   

    def readfile(self, filename):
        if not os.path.exists(filename):
            common.terminate("No such file %s" % filename)

        data = loadtxt(filename , dtype = str, delimiter=self.DELIMITER)

        if data.ndim != 2:
            common.terminate("Something wrong with the data in file %s. It is not 2D matrix" % filename)

        return data

    def get_value_by_title(self, data, data_titles, query_titles, data_type):
        indices = []

        # few queries titles are asked 
        if type(query_titles) == list:
            for query_title in query_titles:
                index_val = where(data_titles == query_title)[0] 
                if index_val:
                    indices.append(index_val[0])
            if indices:
                return data[:, indices].astype(data_type)
        
        # only one query title was asked for
        else:
            query_title = query_titles
            index_val = where(data_titles == query_title)[0] 
            if index_val:
                return data[:, index_val[0]].astype(data_type)
        
        # no title was found
        return None

    def parsedata(self, data):
        titles = data[0,:]
        test_name = titles[0].split(':')[0]
        data = data[1:, :]
        values_number = data.shape[1]
        # obligatory values
        cpgnames = data[:, self.CPGNAMES_INDEX]
        chromosomes = data[:, self.CHR_INDEX]
        positions = data[:, self.POSITION_INDEX]
        pvalues = data[:, self.PVALUE_INDEX].astype(float)
        qvalues = data[:, self.QVALUE_INDEX].astype(float)
        genes = data[:, self.GENE_INDEX]
        categories = data[:, self.CATEGORY_INDEX]

        # optional values
        intercept =        self.get_value_by_title(data, titles, self.INTERCEPT_TITLE, float)
        covariates_betas = self.get_value_by_title(data, titles, self.get_covars_titles(values_number), float) #number of covariates is maximun the number of values (fields) in the data  file
        beta =             self.get_value_by_title(data, titles, self.BETA_TITLE, float)
        stats =            self.get_value_by_title(data, titles, self.STATS_TITLE, float)
        sigma_e =          self.get_value_by_title(data, titles, self.SIGMA_E_TITLE, float)
        sigma_g =          self.get_value_by_title(data, titles, self.SIGMA_G_TITLE, float)


        return test_name, cpgnames, pvalues, \
               qvalues, stats, intercept, covariates_betas, \
               beta, sigma_g , sigma_e, \
               sitesinfo.SitesInfo(cpgnames, chromosomes, positions, genes, categories)
        
        

