from numpy import dot, sqrt


def low_rank_approximation(A, B, i):
    """
    calculates low rank approximation between the first i columns of two n X m matrixes A and B 
    Note: Both A and B dimensions are n X m 
    """
    return dot(A[:,0:i], B[:,0:i].transpose())

def euclidean_distance(A, B):
    """
    calculates euclidean distance between two n X m matrixes A and B
    Note: Both A and B dimensions are n X m 
    """
    return sqrt(((A - B)**2).sum(axis=0))