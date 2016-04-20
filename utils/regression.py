from sklearn import linear_model
from numpy import column_stack, ones, empty_like
from utils import common

class LinearRegression(object):
    def __init__(self, y, x):
        if y.ndim == 1: # y is a vector
            self.residuals = self._regress_out_from_vector(y, x)
        elif y.ndim == 2: # y is a matrix
            self.residuals = self._regress_out_from_matrix(y, x)

        # TODO Elior, change this
        self.p_value = 0.1
        self.t_statistic = 0.2
        self.coef = 0.3

    def _regress_out_from_vector(self, y, x): 
        """
        Y is an n*1 vector (numpy's ndarray)
        X is an n*d matrix (numpy's ndarray)
        """
        zeros_vector =  ones(len(x))
    
        newx = column_stack((zeros_vector, x))

        regr = linear_model.LinearRegression()
        model = regr.fit(newx, y)
        pred = regr.predict(newx)
        residuals = y - pred
        # coef = regr.coef_
        return residuals

    def _regress_out_from_matrix(self, y, x):
        """
        y and x are matrics
        regress out on y's rows
        """
        regressed = empty_like(y)   
        for i,row in enumerate(y):
            residuals = self._regress_out_from_vector(row, x)
            regressed[i] = residuals

        return regressed

class LogisticRegression(object):
    pass