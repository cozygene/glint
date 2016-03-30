from numpy import dot, sqrt
import pca

def low_rank_approximation(O, k):
    """
    O dimensions are n X m
    """
    pca_out = pca.PCA(O)
    return dot(pca_out.P[:,0:k], pca_out.U[:,0:k].transpose())

def euclidean_distance(A, B):
    """
    calculates euclidean distance between two n X m matrixes A and B
    Note: Both A and B dimensions are n X m 
    """
    return sqrt(((A - B)**2).sum(axis=0))