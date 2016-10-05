import os
import sys
import logging
from sklearn import preprocessing
from numpy import sqrt, savetxt, column_stack, where, delete
from utils import tools, pca, LinearRegression, common
from module import Module

RANKED_FILENAME =       'refactor.rankedlist.txt'
COMPONENTS_FILENAME =   'refactor.components.txt'
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
                  use_covars = None,
                  use_phenos = None,
                  bad_probes_list = [],
                  feature_selection = 'normal',
                  ranked_output_filename = RANKED_FILENAME,
                  components_output_filename = COMPONENTS_FILENAME
                ):
        """
        if use_covars is:
                None - no covariates will take into account
                [] (empty list) - all covariates will take into account
                [covar_name1, ... covar_nameX] - only covar names specified in the list will be takem into account
        """
        self.meth_data = methylation_data
        feature_selection = feature_selection.lower().strip()
        
        # validate and process all variables
        self.use_phenos = use_phenos
        self.use_covars = use_covars
        self.fs_preprocessing =           self._validate_fs(feature_selection) 
        self.k =                          self._validate_k(k)
        self.t =                          self._validate_t(t)
        self.minstd =                     self._validate_stdth(minstd)
        self.num_components =             self._validate_num_comp(num_components)
        self.bad_probes =                 bad_probes_list
        self.ranked_output_filename =     ranked_output_filename
        self.components_output_filename = components_output_filename

    def _validate_fs(self, feature_selection):
        if feature_selection not in self.FEATURE_SELECTION:
            common.terminate("choose feature selection from the options: %s (selected fs: %s)" % ( self.FEATURE_SELECTION, feature_selection ))
        
        if self.use_phenos is None:
            if feature_selection == 'phenotype':
                common.terminate("must provide a phenotype (--pheno) when selected feature 'phenotype'")
            elif feature_selection == 'controls':
                common.terminate("must provide one phenotype (--pheno) in a binary format when selected feature 'controls'")
        elif feature_selection != 'phenotype' and feature_selection != 'controls':
            logging.warning("you selected phenotypes but not a phenotype feature, change selection with the flag --fs (selected fs: %s)" % feature_selection)
            
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

    def run( self ):
        # TODO use module name from config file instead of "ReFACTor"
        logging.info('Starting ReFACTor v%s...' % self.VERSION); 
        self.components, self.ranked_sites_indices = self._refactor()
        self.ranked_sites = self.meth_data.cpgnames[self.ranked_sites_indices]
        logging.info('ReFACTor is Done!')

   

    def _refactor( self ):
        """
        run refactor:
            exclude bad probes
            remove sites with low std
            remove covariates
            run feature selection
            computing the ReFACTor components
            find ranked list of the data features
        """
        self._exclude_bad_probes()
        # self.meth_data.remove_missing_values_sites() # nan are not supported TODO uncomment when supported
        self.meth_data.remove_lowest_std_sites(self.minstd)
        # self.meth_data.replace_missing_values_by_mean() # nan are not supported TODO uncomment when supported
       
        # feature selection
        ranked_list = self._feature_selection()
        logging.info('Computing the ReFACTor components...')
        sites = ranked_list[:self.t] # take best t sites indices
        pca_out = pca.PCA(self.meth_data.data[sites,:].transpose())
        score = pca_out.P

        logging.info('Saving a ranked list of the data features to %s...' % self.ranked_output_filename)
        ranked_list_output = [self.meth_data.cpgnames[index] for index in ranked_list]
        savetxt(self.ranked_output_filename, ranked_list_output, fmt='%s')

        logging.info('Saving the ReFACTor components to %s...' % self.components_output_filename)
        components = score[:,:self.num_components]
        components_output = column_stack((self.meth_data.samples_ids, components))
        savetxt(self.components_output_filename, components_output, fmt='%s')
        
        return components, ranked_list


    def _exclude_bad_probes(self):
        """
        excludes bad sites from data
        """
        logging.info("excluding bad sites...")
        self.meth_data.exclude(self.bad_probes)

    
    def _normal_feature_handler(self, meth_data):
        """
        the normal feature selection does not change the data
        Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
        """
        pass
    
    
    def _phenotype_feature_handler(self, meth_data):
        """
        regress out the phenotype
        Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
        """
        logging.info("Running phenotype feature selection...")
        phenotype = meth_data.get_phenotype_subset(self.use_phenos)
        if phenotype.ndim == 2 and phenotype.shape[1] > 1:
            logging.warning("phenotype feature selection was used with more than one phenotype")
        
        logging.info("regressing out phenotype...")
        meth_data.regress_out(phenotype)
        

    def _controls_feature_handler(self, meth_data):  
        """
        keep only the controls samples in the data
        Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
        """
        logging.info("Running controls feature selection...")
        
        phenotype = meth_data.get_phenotype_subset(self.use_phenos)

        if not tools.is_binary_vector(phenotype):
            common.terminate("phenotype in not a binary vector (must be a 1D binary vector with 'controls' feature selection)")
            
        controls_samples_indices = where(phenotype == 0)[0] # assumes its a 1D binary vector
        if (self.k > len(controls_samples_indices)):
            common.terminate("k cannot be greater than amount of control samples")

        remove_indices = delete(range(meth_data.samples_size), controls_samples_indices)
        meth_data.remove_samples_indices(remove_indices)
        
    def _feature_selection(self):
        fs_meth_data = self.meth_data.copy() # create a copy to work on for feature selection only
        self.fs_preprocessing(fs_meth_data) # 'phenotype' or 'controls' feature selection
        
        covars = fs_meth_data.get_covariates_subset(self.use_covars)
        if covars is not None:
            logging.info("regressing out covariates...")
            fs_meth_data.regress_out(covars)
        else:
            logging.info("ignoring covariates")


        ranked_list = self._calc_low_rank_approx_distances(fs_meth_data) # returns array of the indexes in a sorted order
        del fs_meth_data # no need for the copy of the data after found the ranked list
        return ranked_list

    """
    TODO add explanation
    """
    def _calc_low_rank_approx_distances(self, meth_data): # feature selection
        logging.info('Computing low rank approximation of the input data and ranking sites...')

        x = tools.low_rank_approximation(meth_data.data.transpose(), self.k)
        x = x.transpose()

        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(meth_data.data.transpose()).transform(meth_data.data.transpose())
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x.transpose()).transform(x.transpose())

        # normalization
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )

        # find the distance of each site from its low rank approximation.
        distances = tools.euclidean_distance(An, Bn)

        return  distances.argsort() # returns array of the indexes in a sorted order (the original indexes of the values if the array was sorted) 
