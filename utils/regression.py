from sklearn import linear_model
from sklearn import feature_selection
from numpy import column_stack, ones, zeros, empty_like
from numpy import array
import numpy as np
from numpy.linalg import inv
from scipy.stats import t
import statsmodels.api as sm
import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
# TODO make sure the functions here don't change the values of the data: http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html#sklearn.linear_model.LinearRegression.fit (see the "copy_X" param)

def get_dim(vector):
    if vector.ndim == 1 or (vector.ndim == 2 and vector.shape[1] == 1):
        return 1
    return 2

class LogisticRegression(object):
    def __init__(self):
        pass   

    @staticmethod
    def fit_model(y, x, covars = None):
        """
        y is n X 1 - phenotype
        x is n X 1 - site under test
        covars (optional) is n X p

        Returns three arrays of (1+p+m) X 1 - coefficients, t-statistic and p-values:
                coefficients - the oefficients array where coefficients[0] if the coef of the intercept
                                                           coefficients[-1] if the coef of the site under test (the m from input x)
                                                           coefficients[1],..., coefficients[p+1] the coefficient of the covariates
                t-statistic - again index 0 if for the intercept, index -1 for site under test and 1 to p+1 for covars
                p-values - again index 0 if for the intercept, index -1 for site under test and 1 to p+1 for covars
        to sum up - in order to get the coeffs,  p-values and the t-statistic of the site under test (input x) extract coefficients[-1], t-statistic[-1] and p-values[-1]
        """
        if x.ndim == 1:
            x = x.reshape(-1,1) # make sure dim is (n,1) and not(n,)
        if y.ndim == 1:
            y = y.reshape(-1, 1)

        # X should have a column of ones, the site of interest and the covariates
        X = x
        if covars is not None:
            X = column_stack((covars, X))
        n = X.shape[0] # number of sites
        X = np.concatenate((np.ones((n,1)), X), axis=1)
        
       
        logit = sm.Logit(y,X)
        result = logit.fit(disp=False) # False disable the print of "Optimization terminated successfully" message

        #  from doc - 
        # result.params # The parameters of a fitted model  - same as coef if you print result.summary()
        # result.pvalues # p values
        # result.tvalues # Return the t-statistic for a given parameter estimate.
        return result.params, result.tvalues, result.pvalues #coefficients, t-statistic and p-values
        

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
        if get_dim(y) == 2 :
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

    # # the old fit model - it is tested and working
    # @staticmethod
    # def fit_model(y, x, covars = None):
    #     """
    #     y is n X 1, x is n X m and covars (optional) is n X p
    #     Returns three values, each is an m length vector, one for the coefficients, one for the f-statistics and one for the p-values.
    #     """
    #     # Note:  regr.fit(x,y) expects x of 2dim and y of 2 dim
    #     #        f_regression(x,y) expects x of 2dim and y of 1 dim
    #     # therefor we create a 2dim and ydim vectors of x and y and we call them:
    #     # x_res2 - for 2 dim x
    #     # x_res1 - for 1dim x (if possible)
    #     # and same for y
        

    #     # regress out covariates
    #     if covars is not None:
    #         y = LinearRegression.regress_out(y, covars)

    #     # create a 2dim (n,1) and 1dim (n,) vectors of y since there are function expects 2dim and functions expects 1dim
    #     if y.ndim == 1:
    #         y_res1 = y
    #         y_res2 = y.reshape(-1, 1)
    #     else:
    #         y_res1 = y.reshape(-1)
    #         y_res2 = y
        
    #     regr = linear_model.LinearRegression(True)

        
    #     if get_dim(x) == 1: # if x is n X 1
    #         m = 1
    #         if covars is not None:
    #             x_res = LinearRegression.regress_out(x, covars)
    #         else:
    #             x_res = x

    #         # create a 2dim (n,1) and 1dim (n,) vectors of x
    #         if x_res.ndim == 1:
    #             x_res2 = x_res.reshape(-1,1)
    #         else:
    #             x_res2 = x_res

    #         regr.fit(x_res2, y_res2)
    #         coefs = regr.coef_[0]
    #         fstats, pvals = feature_selection.f_regression(x_res2, y_res1, center=True)

    #     else: 
    #         m = x.shape[1]
    #         pvals = zeros(m)
    #         fstats = zeros(m)
    #         coefs = zeros(m)
    #         for i in range(m):
    #             x_res2 = x[:,i].reshape(-1, 1)
    #             if covars is not None:
    #                 x_res1 = LinearRegression.regress_out(x_res2, covars) # regress_out of (n,1) or (n,) vector returns (n,) vector
    #                 x_res2 = x_res1.reshape(-1,1)
    #             regr.fit(x_res2, y_res2)
    #             coefs[i] = regr.coef_[0]
    #             fstats[i], pvals[i] = feature_selection.f_regression(x_res2, y_res1, center=True)

    #     return coefs, fstats, pvals


    @staticmethod
    def fit_model(y, x, covars = None):
        """
        y is n X 1 - phenotype
        x is n X 1 - site under test
        covars (optional) is n X p

        Returns three arrays of (1+p+m) X 1 - coefficients, t-statistic and p-values:
                the first is the coefficients array where coefficients[0] if the coef of the intercept
                                                          coefficients[-1] if the coef of the site under test (the m from input x)
                                                          coefficients[1],..., coefficients[p+1] the coefficient of the covariates
                the second array holds the f-statistics - again index 0 if for the intercept, index -1 for site under test and 1 to p+1 for covars
                the third array holds the p-values - again index 0 if for the intercept, index -1 for site under test and 1 to p+1 for covars
        to sum up - in order to get thecoeffs,  p-values and the t-statistic of the site under test (input x) extract coefficients[-1], t-statistic[-1] and p-values[-1]
        """
        if x.ndim == 1:
            x = x.reshape(-1,1) # make sure dim is (n,1) and not(n,)
        if y.ndim == 1:
            y = y.reshape(-1, 1)

        X = x
        if covars is not None:
            X = column_stack((covars, X))
        
        regr = linear_model.LinearRegression(False)
        n = X.shape[0] # number of sites
        X = np.concatenate((np.ones((n,1)), X), axis=1)

        mdl = regr.fit(X,y)
        sse = np.sum((mdl.predict(X) - y) ** 2, axis=0) / float(X.shape[0] - X.shape[1])
        se = np.array([
            np.sqrt(np.diagonal(sse[i] * np.linalg.inv(np.dot(X.T, X))))
                                                    for i in range(sse.shape[0])
                    ])

        Ts = mdl.coef_ / se
        p = 2 * (1 - t.cdf(np.abs(Ts), y.shape[0] - X.shape[1]))
        return mdl.coef_.reshape(-1), Ts.reshape(-1), p.reshape(-1)  #coefficients, t-statistic and p-values

        # if x.ndim == 1:
        #     x = x.reshape(-1,1) # make sure dim is (n,1) and not(n,)
        # if y.ndim == 1:
        #     y = y.reshape(-1, 1)

        # X = x
        # if covars is not None:
        #     X = column_stack((covars, X))
        
        # regr = linear_model.LinearRegression(False)
        # n = X.shape[0] # number of sites
        # X = np.concatenate((np.ones((n,1)), X), axis=1)
        # p = X.shape[1] # number of covars
        # mdl = regr.fit(X,y)
        # beta = mdl.coef_.reshape(-1)# Beta contains the coefficients of the intercept (beta[0]) and the other features

        # C = (np.sum([(y[i]-np.dot(X[i,:],beta))**2 for i in range(n)]) / (n-5)) * inv(np.dot(X.T,X)) #why is there a (n-5)??
        # pvals = np.empty((p,1)).reshape((p,))
        # Ts = np.empty((p,1)).reshape((p,))

        # for i in range(p):
        #     Ts[i] = beta[i]/((C[i,i])**0.5) # The t-statistic
           
        #     pvals[i] = t.sf(abs(Ts[i]), df=n-p)*2
        
        # return beta, Ts, pvals #coefficients, t-statistic and p-values
