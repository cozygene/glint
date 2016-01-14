# Y is an n*1 vector (numpy's ndarray)
# X is an n*d matrix (numpy's ndarray)

from sklearn import linear_model
from numpy import column_stack, ones    

class LinearRegression(object):
    def __init__(self, y, x): # y is vector x matrix
        zeros_vector =  ones(len(x))
    
        newx = column_stack((zeros_vector, x))

        regr = linear_model.LinearRegression()
        model = regr.fit(newx, y)
        pred = regr.predict(newx)
        self.residuals = y - pred
        self.coef = regr.coef_

class LogisticRegression(object):
    pass