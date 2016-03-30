from numpy import corrcoef
def correlation(x, y):
    assert len(x.shape) == 1, "x must be a 1d vector"
    assert len(y.shape) == 1, "y must be a 1d vector"
    cor = corrcoef(x,y)
    return  (1- abs(cor[0][1])) < 1e-4