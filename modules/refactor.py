import logging
from sklearn import preprocessing
from sklearn.decomposition import PCA
from numpy import dot, linalg, sqrt

VERSION = 1.0 #TODO move to other place

def refactor(O, K, t, num_components, ranked_filename='refactor.out.rankedlist.txt', components_filename='refactor.out.components.txt'):

    logging.info('Starting ReFACTor v%s...' % VERSION);

    
    print O
    print O.transpose()
    # sample_id <- O[1, -1] # extract samples ID
    # O <- O[-1,] # remove sample ID from matrix
    # cpgnames <- O[, 1] ## set rownames
    # O <- O[, -1] 
    # O = matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O))

    # num_components = K;

    logging.info('Running a standard PCA...')
    scaler = preprocessing.StandardScaler().fit(O.transpose())
    pca = PCA(n_components=num_components)
    pcs = pca.fit(O.transpose())
    # pcs = prcomp(scale( t(O) ));



    coeff = pcs.components_
    # coeff = pcs$rotation #U   pcs.components_
    score = pcs.transform(O.transpose())
    # score = pcs$x #P pcs.transform(X)

    logging.info('Compute a low rank approximation of input data and rank sites...')
    x = dot(score[:,1:K], coeff[:,1:K].transpose())
    # x = score[,1:K]%*%t(coeff[,1:K]);
    An = preprocessing.StandardScaler().fit( O.transpose(), with_mean = True, with_std = False )
    # An = scale(t(O),center=T,scale=F)
    Bn = preprocessing.StandardScaler().fit( x, with_mean = True, with_std = False )
    # Bn = scale(x,center=T,scale=F)

    An = ( An.transpose() * (1 / sqrt(numpy.apply_along_axis(sum, 0, linalg.matrix_power(An,2))))).transpose() 
    # An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))

    Bn = ( Bn.transpose() * (1 / sqrt(numpy.apply_along_axis(sum, 0, linalg.matrix_power(Bn,2))))).transpose() 
    # Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))


    # # Find the distance of each site from its low rank approximation.
    distances = numpy.apply_along_axis(sum, 0, linalg.matrix_power(An - Bn,2))
    # distances = apply((An-Bn)^2,2,sum)^0.5 ;

    dsort = sort(distances) ???
    # dsort = sort(distances,index.return=T);
    ??? https://stat.ethz.ch/R-manual/R-devel/library/base/html/sort.html
    # ranked_list = dsort$ix

    logging.info('Compute ReFACTor components...')
    sites = ranked_list[1:t]
    # sites = ranked_list[1:t];
    scaler = preprocessing.StandardScaler().fit(O[sites,:].transpose())
    pcs = PCA(n_components=num_components).fit(O[sites,:].transpose())
    pcs = prcomp(scale(t(O[sites,])));
    # pcs = prcomp(scale(t(O[sites,])));
    first_score = score[:,1:K]
    # first_score <- score[,1:K];
    score = pcs.transform(O.transpose())
    # score = pcs$x

    logging.info('Saving a ranked list of the data features...');
    file(ranked_filename, 'wb').write( ??? )
    # write(t(cbind(ranked_list,cpgnames[ranked_list])),file=ranked_filename,ncol=2)

    logging.info('Saving the ReFACTor components...');
    file(components_filename, 'wb').write( ??? )
    # write(t(score[,1:num_components]),file=components_filename,ncol=num_components)

    return score[,1:num_components], ranked_list, first_score 
    # result <- list(refactor_components=score[,1:num_components],ranked_list=ranked_list, first_pca = first_score) 
    # return(result)


