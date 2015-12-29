import os
import sys
import logging
from sklearn import preprocessing
from numpy import dot, linalg, sqrt, hstack, loadtxt, empty_like
from utils import tools, pca, LinearRegression
from module import Module

#TODO NOTE remember to copy the matrix before making changes!!!!

RANKED_FILENAME =       'refactor.out.rankedlist.txt'
COMPONENTS_FILENAME =   'refactor.out.components.txt'

class Refactor( Module ):
    VERSION = 1.0 #TODO move to other place?

    # all feature selection options. Note: if you add more options you need to write a handler (function that is named by the FEATURE_FUNC_NAME_FORMAT) for each option and 
    FEATURE_SELECTION = ['normal', 'phenotype', 'controls']
    FEATURE_FUNC_NAME_FORMAT = "_{feature_option_name}_feature_handler"   # feature selections function name format

    def __init__( self,
                  methylation_data,
                  K,
                  t = 500,
                  num_components = None, 
                  phenofile = None,
                  covar = None,
                  feature_selection = 'normal',
                  ranked_output_filename = RANKED_FILENAME,
                  components_output_filename = COMPONENTS_FILENAME
                ):

        """
        methylation_data is a MethylationData object
        TODO add class doc here
        """
        self.meth_data = methylation_data
        self.meth_data.data = self.meth_data.data
        feature_selection = feature_selection.lower().strip()
        
        # validate and process all variables
        self.phenotype =                  self._validate_phenotype(phenofile, feature_selection)
        self.covar =                      self._validate_covar(covar)
        self.feature_selection_handler =  self._validate_fs(feature_selection) 
        self.k =                          self._validate_k(K)
        self.t =                          self._validate_t(t)
        self.num_components =             self._validate_num_comp(num_components)
        self.ranked_output_filename =     self._validate_filepath(ranked_output_filename)
        self.components_output_filename = self._validate_filepath(components_output_filename)

    def _validate_phenotype(self, phenofile, feature_selection):
        pheno = None
        if phenofile:
            pheno = self._validate_matrix_ids_and_reorder(phenofile)
            if len(pheno[0]) != 2:
                print("the phenotype file provided is not in the right format. should be 2 columns: 1 - sample id, 2 - phenotype") #TODO is this right?
                sys.exit(2)

            pheno = pheno[:,1:].astype(float) # TODO should check if can convert  to float
            
            if feature_selection == 'normal':
                print("[WARNING] provided phenotype file without fs") #TODO should we warn?

        return pheno

    def _validate_covar(self, covariates):
        covar = None
        if not covariates:
            print("[WARNING] didn't supply covariates file")
        else:
            covar = self._validate_matrix_ids_and_reorder(covariates)
            if len(covar[0]) < 2:
                print("the covariates file provided is not in the right format. should be at least 2 columns") #TODO is this right?
                sys.exit(2)

            covar = covar[:,1:].astype(float)

        return covar

    """
    reads matrix from matrix_file_path
    validates that the matrix has number of rows as the number of sample ids
    checks that the sample ids in matrix (the first column) are the same ids as in sample_ids list
    if they are the same but in different order, reorder matrix rows as in sample_ids
    """
    def _validate_matrix_ids_and_reorder(self, matrix_file_path):
        if not os.path.exists(matrix_file_path) :
            print("The file '%s' doesn't exist. Exiting" % matrix_file_path)
            sys.exit(2)

        data = loadtxt(matrix_file_path, dtype = str)
        if len(data) != len(self.meth_data.samples_ids):
            print("the file provided %s doesn't include all sample ids" % matrix_file_path)
            sys.exit(2)

        matrix_sample_ids = data[:,0]

        if not (self.meth_data.samples_ids == matrix_sample_ids).all(): #todo check this is not by order
            if len(set(self.meth_data.samples_ids)^set(matrix_sample_ids)) != 0:
                print("sample ids in phenotype file are not the same as in the data file")
                sys.exit(2)
            
            print("sample ids in phenotype file are not in the same order as in the datafile, reordering...")
            sample_to_row = dict()
            for sample in data:
                sample_to_row[sample[0]] = sample

            orderd_data = empty_like(data)  
            for i,sid in enumerate(self.meth_data.samples_ids):
                orderd_data[i,:] = sample_to_row[sid]
            return orderd_data

        return data

    """
    must be called after _validate_phenotype
    """
    def _validate_fs(self, feature_selection):
        if feature_selection not in self.FEATURE_SELECTION:
            print("choose fs from feature_selection options: %s (selected fs: %s)" % ( self.FEATURE_SELECTION, feature_selection ))
            sys.exit(2)
        elif feature_selection != 'normal' and self.phenotype is None:
            print("must provide a phenotype file when selected feature 'controls'")
            sys.exit(2)
        elif feature_selection == 'controls' and not self._is_binary_vector(self.phenotype):
            print("must provide a phenotype file in a binary format when selected feature 'controls'")
            sys.exit(2)

        return getattr(self, self.FEATURE_FUNC_NAME_FORMAT.format(feature_option_name=feature_selection))

    """
    2 <= K <= samples size
    """
    def _validate_k(self,k):
        if not (k >= 2 and k <= self.meth_data.samples_size):
            print("k must be at least 2 and smaller than samples size. k = %s, samples = %s" % (k, self.meth_data.samples_size))
            sys.exit(2) 

        return k

    """
    K <= t <= sites size
    must be called after _validate_k
    """
    def _validate_t(self,t):
        if t > self.meth_data.sites_size or t < self.k : 
            print("t cannot be greater than number of sites or smaller than k . t = %s, sites = %s, k = %s" % (t, self.meth_data.sites_size, self.k))
            sys.exit(2) 

        return t

    """
    K <= num_comp  <= samples size
    must be called after _validate_k
    """
    def _validate_num_comp(self,num_comp):
        if num_comp and not (num_comp >= self.k and num_comp <= self.meth_data.samples_size):
            print("number of components must be at least k and smaller than samples size. num_comp = %s, samples = %s, k = %s" % (t, self.meth_data.samples_size, self.k))
            sys.exit(2) 

        return num_comp if num_comp else self.k

    def _validate_filepath(self, filepath):
        # TODO add a varification, try to open the file? 
        return filepath

    """
    gets a vector of ints/doubles and checks:
        if that is a vector
        if all it's values are 0 or 1
    """
    def _is_binary_vector(self, vector):
        values = self.phenotype.squeeze()
        if len(values.shape) != 1: # two dimentions is not a vector
            return False

        if set(values) != set([0,1]): #has values that are not 0 or 1
            return False       

        return True

    def run( self ):
        logging.info('Starting ReFACTor v%s...' % self.VERSION);
        self.components, self.ranked_sites = self._refactor()

   
    """
    writes data to file filepath.
    removes the file if already exists
    """
    def _write_file( self, filepath, data):   
        if  os.path.exists(filepath):
            os.remove(filepath) #TODO dont remove.. do other thing..!
        with open(filepath, 'w') as f:
            f.write(data)

    """
    TODO add doc
    """
    def _refactor( self ):
        self.meth_data.remove_missing_values_sites()
        self.meth_data.remove_lowest_variance_sites()
        self.meth_data.replace_missing_values_by_mean()
        self._remove_covariates()

        # feature selection
        O_find_best_sites = self.feature_selection_handler() #TODO change name of O_find_best_sites
        ranked_list = self._calc_low_rank_approx_distances(O_find_best_sites, self.k)

        logging.info('Compute ReFACTor components...')
        sites = ranked_list[0:self.t]

        pca_out = pca.PCA(self.meth_data.data[sites,:])
        score = pca_out.P

        logging.info('Saving a ranked list of the data features...')
        data = '\n'.join(['%s\t%s'% (index, self.meth_data.cpgnames[index]) for index in ranked_list])
        self._write_file(self.ranked_output_filename, data)

        logging.info('Saving the ReFACTor components...')
        data = '\n'.join(['\t'.join([str(i) for i in line]) for line in score[:,0:self.k]])
        self._write_file(self.components_output_filename, data)
        
        return score[:,0:self.num_components], ranked_list

    def _remove_covariates(self):
        if self.covar is not None:
            O_reg = empty_like(self.meth_data.data)   
            for i,site in enumerate(self.meth_data.data):
                lin_reg = LinearRegression(site, self.covar)
                if len(lin_reg.residuals) != len(site):
                    print "ERROR" #TODO move to test?
                O_reg[i] = lin_reg.residuals
            self.meth_data.data = O_reg


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
        O_tag = empty_like(self.meth_data.data)
        for i,site in enumerate(self.meth_data.data):
            lin_reg = LinearRegression(site, self.phenotype)
            if len(lin_reg.residuals) != len(site): #TODO move to test?
                print "ERROR"
            O_tag[i] = lin_reg.residuals

        return O_tag
    
    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _controls_feature_handler( self ):   
        controls_samples_indices = [i for i, control in enumerate(self.phenotype) if control == 0]
        if (self.k > controls_samples_indices):
            print("k cannot be greater than controls sample")
            sys.exit(2)
        self.meth_data.data = self.meth_data.data[:, controls_samples_indices]
        return self.meth_data.data


    """
    TODO add explanation
    """
    def _calc_low_rank_approx_distances( self, O, i ):

        logging.info('Running a standard PCA...')
        # PCA calc is here so we don't need to pass P and U as arguments to this function. That way this func can be called from feature selection    pca_out = pca.PCA(self.meth_data.data.transpose()) 
        pca_out = pca.PCA(O.transpose()) 

        logging.info('Compute a low rank approximation of input data and rank sites...')
        x = tools.low_rank_approximation(pca_out.P, pca_out.U, i)
       
        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(O.transpose()).transform(O.transpose()) #TODO move transpose out?
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x).transform(x)
        # normalization
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )

        # find the distance of each site from its low rank approximation.
        distances = tools.euclidean_distance(An, Bn)

        return  distances.argsort() # returns array of the indexes in a sorted order (the original indexes of the values if the array was sorted) 
