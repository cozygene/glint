import os
import logging
from sklearn import preprocessing
from sklearn.decomposition import PCA
from numpy import dot, linalg, sqrt,hstack

VERSION = 1.0 #TODO move to other place

def refactor(O, K, t, num_components, ranked_filename='refactor.out.rankedlist.txt', components_filename='refactor.out.components.txt'):

    logging.info('Starting ReFACTor v%s...' % VERSION);
  
    sample_ids = O[0,:][1:] # extract samples ID
    cpgnames = O[:,0][1:]   # extract methylation sites names
    O = O[1:,1:].astype(float) # remove sample ID and cpgnames from matrix
    
    logging.info('Running a standard PCA...')
    scaledO = preprocessing.StandardScaler().fit(O.transpose()).transform(O.transpose())
    pca = PCA().fit(scaledO) 

    coeff = pca.components_.transpose()
    score = pca.transform( scaledO )
    
    logging.info('Compute a low rank approximation of input data and rank sites...')
    x = dot(score[:,0:K] , coeff[:,0:K].transpose())
    An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit( O.transpose() ).transform(O.transpose())
    Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit( x ).transform(x)
    An = (An * 1 / sqrt((An**2).sum(axis=0)))
    Bn = ( Bn * (1 / sqrt((Bn**2).sum(axis=0))))

    # Find the distance of each site from its low rank approximation.
    distances = sqrt(((An - Bn)**2).sum(axis=0))#**0.5
    ranked_list = distances.argsort()
    distances.sort() 

    logging.info('Compute ReFACTor components...')
    sites = ranked_list[0:t]
    scaledOsites = preprocessing.StandardScaler().fit(O[sites,:].transpose()).transform(O[sites,:].transpose())
    pca = PCA().fit(scaledOsites)
    first_score = score[:,0:K]
    score = pca.transform(scaledOsites)

    logging.info('Saving a ranked list of the data features...');

    if  os.path.exists(ranked_filename):
        os.remove(ranked_filename) #TODO  dont remove.. do other thing..!
    with open(ranked_filename, "a") as f:
        for index in ranked_list:
            f.write("%s %s\n"% (index, cpgnames[index]))

    logging.info('Saving the ReFACTor components...');

    if  os.path.exists(components_filename):
        os.remove(components_filename) #TODO dont remove.. do other thing..!
    with open(components_filename, "a") as f:
        for line in score[:,0:K]:
            f.write("\t".join([str(i) for i in line]) + "\n")
    
    return score[:,0:num_components], ranked_list, first_score 

