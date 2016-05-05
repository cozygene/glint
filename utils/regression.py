from sklearn import linear_model
from sklearn import feature_selection
from numpy import column_stack, ones, zeros, empty_like
from numpy import array

# TODO make sure the functions here don't change the values of the data: http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html#sklearn.linear_model.LinearRegression.fit (see the "copy_X" param)

def get_dim(vector):
    if vector.ndim == 1 or (vector.ndim == 2 and vector.shape[1] == 1):
        return 1
    return 2

class LinearRegression(object):

    def __init__(self):
        pass       
        
    @staticmethod
    def regress_out(y, x):
        """
        y is n X d, x is n X p
        returns y_res, an n X d matrix such that the i-th column of y_res contains the residuals of the i-th column of y after regressing out x.
        """
        regr = linear_model.LinearRegression(True)
        n = y.shape[0]

        # if x is a n X 1 vector in the format (n,), change it's size to (n,1)
        if x.ndim == 1:
            x = x.reshape(-1,1)
        
        if get_dim(y) ==2 :
            d = y.shape[1]
            y_res = zeros([n,d])
            for i in range(d):
                regr.fit(x, y[:,i].ravel())
                pred = regr.predict(x)
                y_res[:,i] = y[:,i] - pred

        else:
            if y.ndim == 2:
                y = y.reshape(-1,)  # make a (n,1) vector to (n,) vector
            regr.fit(x, y)
            pred = regr.predict(x)
            y_res = y - pred

        return y_res

    @staticmethod
    def fit_model(y, x, covars = None):
        """
        y is n X 1, x is n X m and covars (optional) is n X p
        Returns three values, each is an m length vector, one for the coefficients, one for the f-statistics and one for the p-values.
        """       
        y_res = y# TODO .ravel()
        regr = linear_model.LinearRegression(True)
        if covars is not None:
            y_res = LinearRegression.regress_out(y_res,covars)
        
        if x.ndim == 1 or (x.ndim == 2 and x.shape[1] == 1): # 1 dim
            m = 1
            if x.ndim == 1:
                x_res = x.reshape(-1, 1) # make a (n,) vector to (n,1) vector
            if covars is not None:
                x_res = LinearRegression.regress_out(x_res,covars)
            regr.fit(x_res, y_res)
            coefs = regr.coef_[0]
            fstats, pvals = feature_selection.f_regression(x_res, y_res, center=True)

        else:
            m = x.shape[1]
            pvals = zeros(m)
            fstats = zeros(m)
            coefs = zeros(m)
            for i in range(m):
                x_res = x[:,i].reshape(-1, 1)
                if covars is not None:
                    x_res = LinearRegression.regress_out(x_res,covars)
                regr.fit(x_res, y_res)
                coefs[i] = regr.coef_[0]
                fstats[i], pvals[i] = feature_selection.f_regression(x_res, y_res, center=True)

        return coefs, fstats, pvals




class LogisticRegression(object):
    pass