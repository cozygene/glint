import logging
from numpy import dot, sqrt, diag, argsort, where, log10, inf
from numpy import min as npnim
import pca
from scipy.linalg import eigh
from scipy.stats import ranksums
import common
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

def minusLog10(values):
    """
    if the p-value is 0 then don't apply -log on it and assign a high number (say 50).
    changes values in place and returns it
    """
    values = -log10(values)
    zeros_pvalues_indices = where(values == inf)[0]
    values[zeros_pvalues_indices] = 50
    return values

def low_rank_approximation(O, k):
    """
    O dimensions are n X m
    """
    import time
    a = time.time()
    pca_out = pca.PCA(O)
    b = time.time()
    res = dot(pca_out.P[:,0:k], pca_out.U[:,0:k].transpose())
    c = time.time()
    logging.debug("PCA TOOK %s SECONDS AND DOT(MULTI) TOOK %s SECONDS" % (b-a,c-b))
    return res

def euclidean_distance(A, B):
    """
    calculates euclidean distance between two n X m matrixes A and B
    Note: Both A and B dimensions are n X m 
    """
    return sqrt(((A - B)**2).sum(axis=0)) # another option is " linalg.norm(A-B, axis=0)"



def symmetrize(X):
    """
    creates a symmetric matrix out of a matrix in which the under (or above) diagonal elements are all zero.
    for example:
    input is 
    1 2 3
    0 4 5
    0 0 6 

    output is
    1 2 3
    2 4 5
    3 5 6   
    """
    return X + X.T - diag(X.diagonal())

def eigenDecompose(X):
    """
    X is 2D Hermitian or symmetric matrix
    returns  eigenvalues and eigenvectors of X 
    """
    logging.info("computing eigendecomposition...")
    s,U = eigh(X) #Returns the eigenvalues and eigenvectors
    if (npnim(s) < -1e-4):
        common.terminate("Negative eigenvalues found")
    s[s<0]=0    
    ind = argsort(s)
    ind = ind[s>1e-12]
    U = U[:, ind]
    s = s[ind]
    return s,U


def standardize(X, axis=0):
    """
    returns normalized matrix
    parameters:
    X - is a 2D matrix of dimensions a by b 
    axis - axes along which the means are computed. default is 0 (mean of b)
    
    return the matrix normalized by axis  
    """
    sites_mean = X.mean(axis=axis)
    sites_std = X.std(axis=axis)
    X -= sites_mean
    X /= sites_std
    return X

def FDR(pvalues):
    """
    pvalues - list of p p-values
    returns a list of q-values
    """
    # res[1] is a vector with the "q-values" (these are the FDR-adjusted p-values)
    res = fdrcorrection0(pvalues) 
    return res[1]


def wilcoxon_test(y, x):
    """
    y - a binary vector (phenotype)
    x - site under test

    returns U statistic and p-value
    """
    indices_0 = where(y==0)[0] # not "sick"
    indices_1 = where(y==1)[0] # "sick"
    x0 = x[indices_0]
    x1 = x[indices_1]
    # Calculate the Z statistic and the p-value of sick elemente vs non-sick elements
    zstat ,pval = ranksums(x0, x1)
    return zstat, pval

def is_binary_vector(vector):
    """
    gets a vector of ints/doubles and returns True if:
        if that is a vector and
        if all it's values are 0 or 1
    otherwise returns false
    """
    if not (vector.ndim == 1 or (vector.ndim == 2 and vector.shape[1] == 1)): # two dimentions is not a vector
        return False

    values = vector.squeeze()
    if set(values) != set([0,1]): #has values that are not 0 or 1
        return False       

    return True