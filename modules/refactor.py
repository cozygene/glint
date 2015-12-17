import os
import logging
from sklearn import preprocessing
from numpy import dot, linalg, sqrt,hstack
from utils import tools, pca
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
                  O,
                  K,
                  t = 500,
                  num_components = None, 
                  phenofile = None,
                  covar = None,
                  feature_selection = 'normal',
                  ranked_output_filename = RANKED_FILENAME,
                  components_output_filename = COMPONENTS_FILENAME
                ):

        feature_selection = feature_selection.lower().strip()
        self.validate(O, K, t, num_components, phenofile, covar, feature_selection, ranked_output_filename, components_output_filename)

        self.k = K
        self.t = t
        self.num_components = num_components if num_components else self.k
        self.ranked_output_filename = ranked_output_filename
        self.components_output_filename = components_output_filename

        self.sample_ids = O[0,:][1:] # extract samples ID
        self.cpgnames = O[:,0][1:]   # extract methylation sites names
        self.O = O[1:,1:].astype(float) # remove sample ID and cpgnames from matrix

        # find the right feature selection function that handels this feature_selection option
        self.feature_selection_handler = getattr(self, self.FEATURE_FUNC_NAME_FORMAT.format(feature_option_name=feature_selection))

    def validate( self, O, K, t, num_components, phenofile, covar, feature_selection, ranked_output_filename, components_output_filename ):
        if phenofile and not os.path.exists(phenofile) :
            print("The file '%s' doesn't exist. Exiting" % args.pheno)
            sys.exit(2)

        if feature_selection not in self.FEATURE_SELECTION:
            print("choose - normal (defualt) (the standard refactor feature selection; default value), controls (by using the controls only; possible only if a binary phenotype is provided), phenotype (by using a subspace orthogonal to the phenotype)")
            sys.exit(2)

        if not covar:
            print "WARNING"

        if covar and not os.path.exists(phenofile) :
            print("The file '%s' doesn't exist. Exiting" % covar)
            sys.exit(2) 

        return True

    def run( self ):
        logging.info('Starting ReFACTor v%s...' % self.VERSION);
        self.components, self.ranked_sites = self._refactor()

        # if covar: run linreg and get residuals

   
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
        # feature selection
        ranked_list = self.feature_selection_handler()

        logging.info('Compute ReFACTor components...')
        sites = ranked_list[0:self.t]

        pca_out = pca.PCA(self.O[sites,:])
        score = pca_out.P

        logging.info('Saving a ranked list of the data features...')
        data = '\n'.join(['%s\t%s'% (index, self.cpgnames[index]) for index in ranked_list])
        self._write_file(self.ranked_output_filename, data)

        logging.info('Saving the ReFACTor components...')
        data = '\n'.join(['\t'.join([str(i) for i in line]) for line in score[:,0:self.k]])
        self._write_file(self.components_output_filename, data)
        
        return score[:,0:self.num_components], ranked_list

    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _normal_feature_handler( self ):
        return self._calc_low_rank_approx_distances(self.O, self.k)
    
    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _phenotype_feature_handler( self ):
        pass
    
    """
    TODO add doc
    Note: function name must be of the format FEATURE_FUNC_NAME_FORMAT
    """
    def _controls_feature_handler( self ):
        #take only people with 0 and extract them from the reqular matrix and give it to _normal_handler 
        pass 

    """
    TODO add explanation
    """
    def _calc_low_rank_approx_distances( self, M, i ):

        logging.info('Running a standard PCA...')
        # PCA calc is here so we don't need to pass P and U as arguments to this function. That way this func can be called from feature selection    pca_out = pca.PCA(self.O.transpose()) 
        pca_out = pca.PCA(M.transpose()) 

        logging.info('Compute a low rank approximation of input data and rank sites...')
        x = tools.low_rank_approximation(pca_out.P, pca_out.U, i)
       
        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(M.transpose()).transform(M.transpose()) #TODO move transpose out?
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x).transform(x)
        # normalization
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )

        # find the distance of each site from its low rank approximation.
        distances = tools.euclidean_distance(An, Bn)

        return  distances.argsort() # returns array of the indexes in a sorted order (the original indexes of the values if the array was sorted) 
