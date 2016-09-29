
EWAS
====

glint allows to perform epigenome-wide association study (EWAS) under several different models, as described bellow. All the different models requires using the general arguments desribed below under "General arguments".


General arguments
^^^^^^^^^^^^^^^^^

**--ewas**

Runs an association test on each site in the data. This will output a results file named results.txt.

.. note:: If --ewas is used without additional arguments glint will use linear regression.

.. note:: glint can produce plots based on the results file. For more details read about the --plot argment.


**--pheno**

Selects a phenotype to use in the association test, according to the phenotype names (the headers provided in the phenotypes file using --phenofile).

For example::

	glint.py --datafile datafile.glint --ewas --linreg --pheno p1

will run EWAS using linear regression model with the phenotype p1.

.. note:: Use the argument --phenofile in order to provide phenotypes.

.. note:: Currently glint can perform EWAS only on one phenotype at a time.



**--covar**

Selects covariates to use in the association ntest, according to the covaraite names (the headers provided in the covariates file using --covarfile).

For example::

	glint.py --datafile datafile.glint --ewas --linreg --pheno p1 --covar c1,c2,c3

will run EWAS using linear regression model with the covaraites c1, c2 and c3.


.. note:: Use the argument --filename in order to provide covaraites.


TODO: mention that --covar with no covariate names will use all the available covaraites.

Linear regression
^^^^^^^^^^^^^^^^^

**--linreg**

Performs EWAS on the data using linear regression model. This is the default model for --ewas.

The output file named results.txt includes a list of the sites, sorted by their association p-value. The output file format includes the following columns: ID (CpG identifier), chromosome (chromosome number of the site), MAPINFO (position of the site in the genome), p-value, q-value, intercept , V1 (coefficient of the first covariate),..., Vn (coefficient of the last covaraite, beta (the coefficient of the site under test), statistic (the test statistic), UCSC_RefGene_Name (name of the gene that is closest to this site), Relation_to_UCSC_CpG_Island (category)

For example::

	glint.py --datafile datafile.glint --ewas --linreg --pheno p1

will run EWAS using linear regression model.




Logistic regression
^^^^^^^^^^^^^^^^^^^

**--logreg**

Performs EWAS on the data using logistic regression model. This option requires a binary phenotype (controls are assumed to be coded as '0' and cases as '1').

The output file named results.txt includes a list of the sites, sorted by their association p-value. The output file format is indentical to the once described under --linreg.

For example::

	glint.py --datafile datafile.glint --ewas --logreg --pheno p1

will run EWAS using logistic regression model.




Wilcoxon rank-sum test
^^^^^^^^^^^^^^^^^^^^^^

**--wilc**

Performs EWAS on the data using the non-parameteric Wilcoxon rank-sum text. This option requires a binary phenotype (controls are assumed to be coded as '0' and cases as '1').

The output file named results.txt includes a list of the sites, sorted by their association p-value. The output file format includes the following columns: ID (CpG identifier), chromosome (chromosome number of the site), MAPINFO (position of the site in the genome), p-value, q-value, statistic (the test statistic), UCSC_RefGene_Name (name of the gene that is closest to this site), Relation_to_UCSC_CpG_Island (category)


For example::

	glint.py --datafile datafile.glint --ewas --wilc --pheno p1

will run EWAS using the Wilcoxon rank-sum test.




Linear mixed model (LMM)
^^^^^^^^^^^^^^^^^^^^^^^^

**--lmm**

Performs EWAS on the data using linear mixed model (LMM). This is an implementation of the FaST-LMM algorithm by Lippert et al. [1]_

The output file named results.txt includes a list of the sites, sorted by their association p-value. The output file includes the following columns:  ID (CpG identifiers), chromosome (chromosome number of the site), MAPINFO (position of the site in the genome), p-value, q-value, intercept , V1 (coefficient of the first covariate),..., Vn (coefficient of the last covaraite, beta (the coefficient of the site under test), statistic (the test statistic), sigma-e (an estimate of sigma_e), sigma-g (an estimate of sigma_g), UCSC_RefGene_Name (name of the gene that is closest to this site), Relation_to_UCSC_CpG_Island (category)


**--kinship**

The kinship matrix for modelling the inter-individual similarity in the data that is required for the LMM. glint allows two options:

- User-supplied kinship - users can suplly a text file with samples by samples kinship matrix (with no row or column headers). 
- ReFACTor - the ReFACTor algorithm can be used for constructing the kinship matrix. If this option is used then ReFACTor is executed for selecting the top informative sites in the data. The kinship matrix is then constructed by calculatign the empirical covariance matrix of the samples based on the selected sites.

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno p1 --kinship kinship.txt

will run EWAS using LMM with the kinship matrix specified in the kinship.txt file. Alternatively:

	glint.py --datafile datafile.glint --ewas --lmm --pheno p1 --kinship refactor --k 6

will use the ReFACTor algorithm for constructing the kinship matrix (where 6 is the number of assumed cell types).


.. note:: If the "refactor" option is used then all of the arguments available with the --refactor arguments are also available here.



**--reml**

Allows to indicate whether rstricted maximum likelihood estimation (REML) or maximum likelihood estimation (ML) should be used. The default value is 1 (REML). Alternatively, ML can be selected usign the value 0.

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno p1 --reml 0

will perform EWAS on the data using LMM with ML estimation.


**--norm**

This argument normalizes the covariates (if supplied) before fitting the LMM.

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno p1 --norm

will perform EWAS on the data using LMM after normalizing the covariates.



**--oneld**

This argument allows to fit the log delta parameter in the Fast-LMM model only once for speed-up (instead for each site separately).

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno p1 --oneld

will perform EWAS on the data using LMM with a single value of log detla.



.. [1] Lippert, Christoph, Jennifer Listgarten, Ying Liu, Carl M. Kadie, Robert I. Davidson, and David Heckerman. "FaST linear mixed models for genome-wide association studies." Nature methods 8, no. 10 (2011): 833-835.
