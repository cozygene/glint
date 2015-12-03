import logging
from sklearn import preprocessing
from sklearn.decomposition import PCA
from numpy import dot, linalg, sqrt,apply_along_axis

VERSION = 1.0 #TODO move to other place

def refactor(O, K, t, num_components, ranked_filename='refactor.out.rankedlist.txt', components_filename='refactor.out.components.txt'):

    print('Starting ReFACTor v%s...' % VERSION);
    # logging.info('Starting ReFACTor v%s...' % VERSION);

    print len(O)
    print len(O[0])    
    # print O[0] # in R print(O[1,])
    # print O.transpose()
    # sample_id <- O[1, -1] # extract samples ID
    # O <- O[-1,] # remove sample ID from matrix
    # cpgnames <- O[, 1] ## set rownames
    # O <- O[, -1] 
    # O = matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O))

    # num_components = K;

    print('Running a standard PCA...')
    # logging.info('Running a standard PCA...')
    scaledO = preprocessing.StandardScaler().fit(O.transpose()).transform(O.transpose())

    # print scaler.transform(O.transpose())
    # print(dir(scaler))
    # import pdb
    # pdb.set_trace()
    pca = PCA(n_components=num_components).fit(scaledO)
    # pcs = pca.fit(O.transpose())
    # pcs = prcomp(scale( t(O) ));



    coeff = pca.components_
    print "coeff"

    print len(coeff)
    print len(coeff[0])
    print coeff[0]
    print "_____________________0"
    # coeff = pcs$rotation #U   pcs.components_
    score = pca.transform( scaledO ) # TODO or score = pca.transform(O.transpose())
    print "score"
    print len(score)
    print len(score[0])
    print score[0]
    # score = pcs$x #P pcs.transform(X)
    # print score

    print "_____________________1"
    logging.info('Compute a low rank approximation of input data and rank sites...')
    x = dot(score[:,1:K], coeff[:,1:K].transpose())
    print len(x)
    print len(x[0])
    print x[0]
    # x = score[,1:K]%*%t(coeff[,1:K]);
    # print x
    print "_____________________2"
    An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit( O.transpose() ).transform(O.transpose())
    # An = scale(t(O),center=T,scale=F)
    # print An.transform(O.transpose())
    
    print "_____________________3"
    Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit( x ).transform(x)
    print len(Bn)
    print len(Bn[0])
    print Bn[0]
    # Bn = scale(x,center=T,scale=F)
    # print Bn.transform(x)

    print "_____________________4"


    An = (An * 1 / sqrt((An**2).sum(axis=0))).transpose() 
    # An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
 
    print "_____________________5"
    Bn = ( Bn * (1 / sqrt((Bn**2).sum(axis=0)))).transpose() 
    # Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))
    print len(Bn)
    print len(Bn[0])
    print Bn[0]
    print "_____________________6"

    # # Find the distance of each site from its low rank approximation.
    distances = apply_along_axis(sum, 1, (An - Bn)**2)
    # distances = apply((An-Bn)^2,2,sum)^0.5 ;
    # print distances
    print "_____________________7"
    dsort = sort(distances) #???
    # dsort = sort(distances,index.return=T);
    # ??? https://stat.ethz.ch/R-manual/R-devel/library/base/html/sort.html
    ranked_list = dsort #???
    # ranked_list = dsort$ix

    print('Compute ReFACTor components...')
    # logging.info('Compute ReFACTor components...')
    sites = ranked_list[1:t]
    # sites = ranked_list[1:t];
    scaledOsites = preprocessing.StandardScaler().fit(O[sites,:].transpose()).transform(O[sites,:].transpose())
    pca = PCA(n_components=num_components).fit(scaledOsites)
    # pcs = prcomp(scale(t(O[sites,])));
    first_score = score[:,1:K]
    # first_score <- score[,1:K];
    score = pca.transform(scaledOsites) # score = pca.transform(O[sites,:].transpose())
    # score = pcs$x

    logging.info('Saving a ranked list of the data features...');
    file(ranked_filename, 'wb').write( ranked_list ) #one line seperate with , ??? 
    # write(t(cbind(ranked_list,cpgnames[ranked_list])),file=ranked_filename,ncol=2)

    logging.info('Saving the ReFACTor components...');
    s=""
    for i in len(score):
        s += score[i][1:num_components]
    file(components_filename, 'wb').write( s) #???
    # write(t(score[,1:num_components]),file=components_filename,ncol=num_components)

    return s, ranked_list, first_score #score[,1:num_components], ranked_list, first_score 
    # result <- list(refactor_components=score[,1:num_components],ranked_list=ranked_list, first_pca = first_score) 
    # return(result)


