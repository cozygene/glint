import os
from utils import LinearRegression, tools#, plot
from numpy import column_stack, ones, savetxt, array, insert, vstack, loadtxt, append
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
            coefs, fstats, p_value = LinearRegression.fit_model(self.meth_data.phenotype, site, covars = self.meth_data.covar) #TODO add test
            output.append([self.meth_data.cpgnames[i], p_value[0], fstats[0], coefs[0]  ])

        
        output.sort(key = lambda x: x[1]) # sort output by p-value
        output = array(output)

        if output_filename:
            qqplot_out = output_filename + '_qqplot' # TODO Elior, change this name (qqplot output file name)?
            logging.info("savings results to %s and qq-plot to %s" % (output_filename, qqplot_out))  
            savetxt(output_filename, output, fmt='%s')
            # plot the p-value
            # qqplot = plot.QQPlot(save_file = qqplot_out)
            # qqplot.draw(output[:,1].astype(float), title = "TODO Elior, CHANGE THIS", xtitle="TODO Elior, change this x", ytitle = "TODO Elior, change this y")

            

            qqplot = plot.ManhattanPlot(save_file = "ManhattanPlot-reut")
            qqplot.draw(  1, self.meth_data.cpgnames, output[:,1].astype(float), title = "TODO Elior, CHANGE THIS", xtitle="TODO Elior, change this x", ytitle = "TODO Elior, change this y")


        return output

class EWASResults(object):
    CPGNAMES_INDEX = 0
    CHR_INDEX = 1
    POSITION_INDEX = 2
    PVALUE_INDEX = 3
    QVALUE_INDEX = 4
    EXTRA_START_INDEX = 5
    EXTRA_END_INDEX = -3
    GENE_INDEX = -2
    CATEGORY_INDEX = -1 #island

    CPGNAMES_TITLE = "{test_name}:ID"
    CHR_TITLE = "chromosome"
    POSITION_TITLE = "MAPINFO" #position
    PVALUE_TITLE = "p-values"
    QVALUE_TITLE = "q-values"
    GENE_TITLE = "UCSC_RefGene_Name" #gene
    CATEGORY_TITLE = "Relation_to_UCSC_CpG_Island" #island


    def __init__(self, test_name, cpgnames, pvalues, extradata, extradata_titles, qvalues = None, sites_info_obj = None):
        self.test_name = test_name
        self.cpgnames = cpgnames
        self.pvalues = pvalues
        self.extradata = extradata
        self.extradata_titles = extradata_titles
        

        if qvalues is None:
            self.qvalues = tools.FDR(pvalues)
        else:
            self.qvalues = qvalues
        
        if sites_info_obj is None:   
            self.sites_info = sitesinfo.SitesInfoGenerator(self.cpgnames)
        else:
            self.sites_info = sites_info_obj

class EWASResultsCreator(EWASResults):
    def __init__(self, test_name, cpgnames, pvalues, extradata = None, extradata_titles = None):
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
        super(EWASResultsCreator, self).__init__(test_name, cpgnames, pvalues, extradata, extradata_titles)
        self.output = self.generate_data()

    def generate_data(self):

        titles = insert(self.extradata_titles, 0, [self.CPGNAMES_TITLE, self.CHR_TITLE, self.POSITION_TITLE, self.PVALUE_TITLE, self.QVALUE_TITLE]) #add test name to title
        append(titles, [self.GENE_TITLE, self.CATEGORY_TITLE])
        data = column_stack((self.cpgnames, self.sites_info.chromosomes, self.sites_info.positions, self.pvalues, self.qvalues, self.extradata, self.sites_info.genes, self.sites_info.categories))
        output = vstack((titles, data))
        return output
        
    def save(self, output_filename):
        logging.info("%s results are saved to file %s" %(self.test_name, output_filename))
        savetxt(output_filename, self.output, fmt = "%s")
        
class EWASResultsParser(EWASResults):
    """docstring for EWASResultsParser"""
    def __init__(self, results_filemame):
        data = self.readfile(results_filemame)
        test_name, cpgnames, pvalues, qvalues, extra, extratitles, sites_info = self.parsedata(data)
        super(EWASResultsParser, self).__init__(test_name, cpgnames, pvalues, extra, extratitles, qvalues = qvalues, sites_info_obj = sites_info)

    def readfile(self, filename):
        if not os.path.exists(filename):
            common.terminate("No such file %s" % filename)

        data = loadtxt(filename , dtype = str)

        if data.ndim != 2:
            common.terminate("Something wrong with the data in file %s. It is not 2D matrix" % filename)

        return data

    def parsedata(self, data):
        titles = data[0,:]
        test_name = titles[0].split(':')[0]

        data = data[1:, :]

        chromosomes = data[:, self.CHR_INDEX]
        positions = data[:, self.POSITION_INDEX]
        genes = data[:, self.GENE_INDEX]
        categories = data[:, self.CATEGORY_INDEX]
        cpgnames = data[:, self.CPGNAMES_INDEX]
        pvalues = data[:, self.PVALUE_INDEX]
        qvalues = data[:, self.QVALUE_INDEX]
        extra = data[:, self.EXTRA_START_INDEX:self.EXTRA_END_INDEX]
        return test_name, cpgnames, pvalues, qvalues, extra, titles[self.EXTRA_START_INDEX:self.EXTRA_END_INDEX], SitesInfo(cpgnames, chromosomes, positions, genes, categories)
        

