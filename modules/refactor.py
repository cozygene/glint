import os
import sys
import logging
from sklearn import preprocessing
from numpy import sqrt, savetxt, column_stack
from utils import tools, pca, LinearRegression, common
from module import Module

RANKED_FILENAME =       'refactor.out.rankedlist.txt'
COMPONENTS_FILENAME =   'refactor.out.components.txt'
#NOTE remember to copy the matrix before making changes!!!!

class Refactor(Module):
    VERSION = 1.0 #TODO move to config file


    # all feature selection options. Note: if you add more options you need to write a handler (function that is named by the FEATURE_FUNC_NAME_FORMAT) for each option and 
    FEATURE_SELECTION = ['normal', 'phenotype', 'controls']
    FEATURE_FUNC_NAME_FORMAT = "_{feature_option_name}_feature_handler"   # feature selections function name format

    def __init__( self,
                  methylation_data,
                  k,
                  t = 500,
                  minstd = 0.02,
                  num_components = None,
                  suppress_covars = False,
                  bad_probes_list = [],
                  feature_selection = 'normal',
                  ranked_output_filename = RANKED_FILENAME,
                  components_output_filename = COMPONENTS_FILENAME
                ):

        self.meth_data = methylation_data
        feature_selection = feature_selection.lower().strip()
        
        # validate and process all variables
        self.feature_selection_handler =  self._validate_fs(feature_selection) 
        self.k =                          self._validate_k(k)
        self.t =                          self._validate_t(t)
        self.minstd =                      self._validate_stdth(minstd)
        self.num_components =             self._validate_num_comp(num_components)
        self.bad_probes =                 bad_probes_list
        self.ranked_output_filename =     ranked_output_filename
        self.components_output_filename = components_output_filename
        self.remove_covariates =          not suppress_covars


    def _validate_fs(self, feature_selection):
        if feature_selection not in self.FEATURE_SELECTION:
            common.terminate("choose fs from feature_selection options: %s (selected fs: %s)" % ( self.FEATURE_SELECTION, feature_selection ))
        elif feature_selection == 'phenotype' and self.meth_data.phenotype is None:
            common.terminate("must provide a phenotype file when selected feature 'phenotype'")
        elif feature_selection == 'controls' and (self.meth_data.phenotype is None or not self._is_binary_vector(self.meth_data.phenotype)):
            common.terminate("must provide a phenotype file in a binary format when selected feature 'controls'")

        return getattr(self, self.FEATURE_FUNC_NAME_FORMAT.format(feature_option_name=feature_selection))

    """
    2 <= k <= samples size
    """
    def _validate_k(self,k):
        if not (k >= 2 and k <= self.meth_data.samples_size):
            common.terminate("k must be at least 2 and smaller than the number of samples size. k = %s, samples = %s" % (k, self.meth_data.samples_size))

        return k

    """
    k <= t <= sites size
    must be called after _validate_k
    """
    def _validate_t(self,t):
        if t > self.meth_data.sites_size or t < self.k : 
            common.terminate("t cannot be greater than the number of sites or smaller than k . t = %s, sites = %s, k = %s" % (t, self.meth_data.sites_size, self.k))

        return t

    """
    0 < stdth < 1
    """
    def _validate_stdth(self, stdth):
        if stdth > 1 or stdth < 0:
            common.terminate("minstd cannot be greater than 1 and smaller than 0. stdth = %s" % stdth)
        return stdth

    """
    k <= num_comp  <= samples size
    must be called after _validate_k
    """
    def _validate_num_comp(self,num_comp):
        if num_comp and not (num_comp >= self.k and num_comp <= self.meth_data.samples_size):
            common.terminate("number of components must be at least k and smaller than the number of samples size. num_comp = %s, samples = %s, k = %s" % (num_comp, self.meth_data.samples_size, self.k))

        return num_comp if num_comp else self.k


    """
    gets a vector of ints/doubles and checks:
        if that is a vector
        if all it's values are 0 or 1
    """
    def _is_binary_vector(self, vector):
        values = self.meth_data.phenotype.squeeze()
        if values.ndim != 1: # two dimentions is not a vector
            return False

        if set(values) != set([0,1]): #has values that are not 0 or 1
            return False       

        return True

    def run( self ):
        # TODO use module name from config file instead of "ReFACTor"
        logging.info('Starting ReFACTor v%s...' % self.VERSION); 
        self.components, self.ranked_sites_indices = self._refactor()
        self.ranked_sites = self.meth_data.cpgnames[self.ranked_sites_indices]
        logging.info('ReFACTor is Done!')

   
    """
    writes data to file filepath.
    removes the file if already exists
    """
    def _write_file( self, filepath, data):   
        if  os.path.exists(filepath):
            os.remove(filepath) #TODO move the last file to other place?
        with open(filepath, 'w') as f:
            f.write(data)

    """
    TODO add doc
    """
    def _refactor( self ):
        self._exclude_bad_probes()
        # self.meth_data.remove_missing_values_sites() # nan are not supported TODO uncomment when supported
        self.meth_data.remove_lowest_std_sites(self.minstd)
        # self.meth_data.replace_missing_values_by_mean() # nan are not supported TODO uncomment when supported
        self._remove_covariates()

        # feature selection
        distances = self.feature_selection_handler()
        ranked_list = self._calc_low_rank_approx_distances()

        logging.info('Computing the ReFACTor components...')
        sites = ranked_list[0:self.t]

        pca_out = pca.PCA(self.meth_data.data[sites,:].transpose())
        score = pca_out.P

        logging.info('Saving a ranked list of the data features...')
        ranked_list_output = [[index + 1, self.meth_data.cpgnames[index]] for index in ranked_list]
        savetxt(self.ranked_output_filename, ranked_list_output, fmt='%s')

        logging.info('Saving the ReFACTor components...')
        components = score[:,0:self.num_components]
        components_output = column_stack((self.meth_data.samples_ids, components))
        savetxt(self.components_output_filename, components_output, fmt='%s')
        
        return components, ranked_list

    def _exclude_bad_probes(self):
        """
        excludes bad sites from data
        """
        logging.info("excluding bad sites...")
        self.meth_data.exclude(self.bad_probes)

    def _remove_covariates(self):
        if self.remove_covariates and self.meth_data.covar is not None:
            logging.info("Removing covariates...")

            lin_reg = LinearRegression(self.meth_data.data, self.meth_data.covar)
            self.meth_data.data = lin_reg.residuals

    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _normal_feature_handler( self ):
        return self.meth_data.data
    
    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _phenotype_feature_handler( self ):
        logging.info("Running phenotype feature selection...")
        
        lin_reg = LinearRegression(self.meth_data.data, self.meth_data.phenotype)
        return lin_reg.residuals

    
    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _controls_feature_handler( self ):  
        logging.info("Running controls feature selection...")
        controls_samples_indices = [i for i, control in enumerate(self.meth_data.phenotype) if control == 0]
        if (self.k > controls_samples_indices):
            common.terminate("k cannot be greater than controls sample")

        self.meth_data.data = self.meth_data.data[:, controls_samples_indices]
        return self.meth_data.data


    """
    TODO add explanation
    """
    def _calc_low_rank_approx_distances(self):
        logging.info('Computing low rank approximation of the input data and ranking sites...')

        x = tools.low_rank_approximation(self.meth_data.data.transpose(), self.k)
        x = x.transpose()

        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(self.meth_data.data.transpose()).transform(self.meth_data.data.transpose())
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x.transpose()).transform(x.transpose())

        # normalization
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )

        # find the distance of each site from its low rank approximation.
        distances = tools.euclidean_distance(An, Bn)

        return  distances.argsort() # returns array of the indexes in a sorted order (the original indexes of the values if the array was sorted) 
