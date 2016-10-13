from numpy import corrcoef
def correlation(x, y, epsilon = 1e-4):
    """
    x,y are ndarrays of dimensions 1
    returns True if x,y correative
    returns False otherwise
    """
    assert x.ndim == 1, "x must be a 1d vector"
    assert y.ndim == 1, "y must be a 1d vector"
    if (x==y).all(): # in case of zeroz arrays
    	return True
    cor = corrcoef(x,y)
    return  (1- abs(cor[0][1])) < epsilon