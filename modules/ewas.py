import os
from utils import LinearRegression, tools#, plot
from numpy import column_stack, ones, savetxt, array, insert, vstack, loadtxt, append, where
from module import Module
from utils import common, plot, sitesinfo
import logging

"""
copy meth_data in advance
"""
class EWAS(Module):
    AVALIABLE_TESTS = ['linear_regression', 'logistic_regression']
    TEST_FUNC_NAME_FORMAT = "_{test_name}_test"   # feature selections function name format

    def __init__(self, methylation_data, tests_list):
        self.meth_data = methylation_data
        self.test_handlers = self._get_test_handler(tests_list)

    def run(self):
        logging.info('starting EWAS...');
        #running association tests
        results = [test_handler(output_filename = 'ewas_' + test_name) for (test_name,test_handler) in self.test_handlers]
        logging.info('EWAS is Done!')
        return results

    def _get_test_handler(self, tests_list):
        # check that the tests in test_list are all optional tests (found in AVALIABLE_TESTS)
        if set(set(tests_list).difference(set(self.AVALIABLE_TESTS))) == 0:
            common.terminate('tests %s are not available' % str(set(tests_list).difference(set(self.AVALIABLE_TESTS))))

        return [(test,getattr(self, self.TEST_FUNC_NAME_FORMAT.format(test_name=test))) for test in tests_list]

    def _logistic_regression_test(self, output_filename = None):
        logging.warning("logistic regression is not supported for the moment...") #todo when implementing remove this
        pass
        
    def _linear_regression_test(self, output_filename = None):
        """
        linear regression test
        """
        logging.info("running linear regression test...")
        output = []           
        for i, site in enumerate(self.meth_data.data):
            coefs, fstats, p_value = LinearRegression.bla(self.meth_data.phenotype, site, covars = self.meth_data.covar) #TODO add test
            # coefs, fstats, p_value = LinearRegression.fit_model(self.meth_data.phenotype, site, covars = self.meth_data.covar) #TODO add test
            #order is: cpgnames, pvalues, t-statistic, intercept coeffs, covariates coeffs (could be a matrix)          
            #order is: cpgnames, pvalues, t-statistic, intercept coeffs, covariates coeffs (could be a matrix), site under test coeffs
            output.append([self.meth_data.cpgnames[i], p_value[-1], fstats[-1], coefs[0], coefs[1:-1], coefs[-1]])

        
        output.sort(key = lambda x: x[1]) # sort output by p-value (1 is p-value index)
        sorted_cpgnames = array(output[:,0])
        sorted_pvalues  = array(output[:,1]).astype(float)
        sorted_fstats   = array(output[:,2]).astype(float)
        sorted_intercept_beta = array(output[:,3]).astype(float)
        sorted_covars_betas   = array(output[:,4]).astype(float)
        sorted_site_beta      = array(output[:,5]).astype(float)
        # output = array(output)

        #cpgnames, pvalues, intercept-beta 
                # intercept_beta,             \
                # covariates_betas,           \
                # site_beta,                  \
                # np.array(sigma_e_est),      \
                # np.array(sigma_g_est),      \
                # np.array(statistics)
        return  sorted_cpgnames , sorted_pvalues, sorted_fstats, sorted_intercept_beta, sorted_covars_betas, sorted_site_beta
        # if output_filename:
            # qqplot_out = output_filename + '_qqplot' # TODO Elior, change this name (qqplot output file name)?
            # logging.info("savings results to %s and qq-plot to %s" % (output_filename, qqplot_out))  
            # savetxt(output_filename, output, fmt='%s')
            # plot the p-value
            # qqplot = plot.QQPlot(save_file = qqplot_out)
            # qqplot.draw(output[:,1].astype(float), title = "TODO Elior, CHANGE THIS", xtitle="TODO Elior, change this x", ytitle = "TODO Elior, change this y")

            

            # qqplot = plot.ManhattanPlot(save_file = "ManhattanPlot-reut")
            # qqplot.draw(  1, self.meth_data.cpgnames, output[:,1].astype(float), title = "TODO Elior, CHANGE THIS", xtitle="TODO Elior, change this x", ytitle = "TODO Elior, change this y")


        return output

    def bla():
        sorted_cpgnames , sorted_pvalues, sorted_fstats, sorted_intercept_beta, sorted_covars_betas, sorted_site_beta = _linear_regression_test()
        num_of_covars = sorted_covars_betas.shape[1]
        covars_beta_titles = ["V%d" % (i+1) for i in range(num_of_covars)]
        additional_results = column_stack((intercept_beta, sorted_covars_betas, site_beta, statistics, sigma_e, sigma_g))
        titles = ['intercept'] + covars_beta_titles + ['beta', 'statistic', 'sigma-e', 'sigma-g']
        
        # generate result - by EWAS output format
        ewas_res = ewas.EWASResultsCreator("LMM", sorted_cpgnames, pvalues, additional_results, array(titles))  

        # save results
        output_file = LMM_OUT_SUFFIX if output_perfix is None else output_perfix + LMM_OUT_SUFFIX
        ewas_res.save(output_file)
        return ewas_res


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
    def __init__(self, test_name, cpgnames, pvalues, qvalues = None, statistic = None, intercept_coefs = None, covars_coefs = None, site_coefs = None, sigma_g = None, sigma_e = None, sites_info_obj = None):
        """
        cpgnames - list of cpgnames
        pvalues - list of pvalues for each cpg (i.e pvalues[i] is the pvalue of site named cpgnames[i])
        extradata - a matrix with extra information about the sites (coefficients, statistics)
                    dimensions of n by t where n is number of sites and t is number of extra data
                    transpose extradata before calling this function if needed
        extradata_titles - the titles to be written in the output file, explaining each column in extradata
                            extradata_title[i] is the title describing the extradata found in extradata[:,i] (i'th column)
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
        
        

