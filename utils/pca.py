from sklearn import preprocessing
from sklearn.decomposition import PCA as pca

class PCA( object ):
    """
    A is centered and each column is devided by its standard deviation before calculating PCA

    input:  A - dimensions are n X m (where m > n)
    output: U - dimensions m X n (loadings)
            P - dimensions n X n (scores)

    Note: A ~ P * U^t
    """
    def __init__(self, A):
        scaled = preprocessing.StandardScaler().fit(A).transform(A)
        pca_res = pca().fit(scaled)
        self.U = pca_res.components_.transpose() # loadings
        self.P = pca_res.transform(scaled)     # scores
