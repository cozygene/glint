


EWAS
====

glint allows to perform epigenome-wide association study (EWAS) under several different models, as described bellow. All of the different models require using the general arguments desribed below under "General arguments".

.. note:: The example commands described bellow assume that the user generated `glint files`_ with covariates file and phenotypes file.


|
|

General arguments
^^^^^^^^^^^^^^^^^

.. _--ewas:

**--ewas**

Runs an association test on each site in the data and outputs a results file.

.. note:: If `--ewas`_ is used without additional arguments glint will use linear regression as a default.

.. note:: glint can produce plots based on the results file. For more details read about the `--plot`_ argument.

.. note:: Add *--out filename* in order to change the default output name.


.. _--pheno:

**--pheno**

Selects a phenotype to use in the association test.

For example::

	glint.py --datafile datafile.glint --ewas --linreg --pheno y1


will run EWAS using linear regression model with the phenotype y1. The names of the phenotypes are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `glint files`_.

.. note:: Use the argument `--phenofile`_ in order to provide phenotypes that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.



**--covar**


Selects covariates to use in the association test.

For example::

	glint.py --datafile datafile.glint --refactor --k 6 --covar c1 c2 c3

will run EWAS using linear regression model with the covariates c1, c2 and c3. The names of the covariates are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `glint files`_.

Alternatively, run::

	glint.py --datafile datafile.glint --refactor --k 6 --covar

without specifying names of covariates in order to include into the model all of the covariates included in the glint file.

.. note:: Use the argument `--covarfile`_ in order to provide covariates that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.


|
|

Linear regression
^^^^^^^^^^^^^^^^^

.. _--linreg:

**--linreg**

Performs EWAS on the data using linear regression model. This is the default model for `--ewas`_.

The output file titled *results.glint.linreg.txt* includes a list of the sites, sorted by their association p-value. The output file format includes the following columns: ID (CpG identifier), chromosome (chromosome number of the site), MAPINFO (position of the site in the genome), p-value, q-value, intercept , V1 (coefficient of the first covariate),..., Vn (coefficient of the last covaraite, beta (the coefficient of the site under test), statistic (the test statistic), UCSC_RefGene_Name (name of the gene that is closest to this site), Relation_to_UCSC_CpG_Island (category)

For example::

	glint.py --datafile datafile.glint --ewas --linreg --pheno y1

will run EWAS using linear regression model.


|
|

Logistic regression
^^^^^^^^^^^^^^^^^^^

**--logreg**

Performs EWAS on the data using logistic regression model. This option requires a binary phenotype (controls are assumed to be coded as '0' and cases as '1').

The output file titled *results.glint.logreg.txt* includes a list of the sites, sorted by their association p-value. The output file format is indentical to the once described under `--linreg`_.

For example::

	glint.py --datafile datafile.glint --ewas --logreg --pheno y1

will run EWAS using logistic regression model.



|
|

Wilcoxon rank-sum test
^^^^^^^^^^^^^^^^^^^^^^

**--wilc**

Performs EWAS on the data using the non-parameteric Wilcoxon rank-sum text. This option requires a binary phenotype (controls are assumed to be coded as '0' and cases as '1').

The output file titled *results.wilc.logreg.txt* includes a list of the sites, sorted by their association p-value. The output file format includes the following columns: ID (CpG identifier), chromosome (chromosome number of the site), MAPINFO (position of the site in the genome), p-value, q-value, statistic (the test statistic), UCSC_RefGene_Name (name of the gene that is closest to this site), Relation_to_UCSC_CpG_Island (category)


For example::

	glint.py --datafile datafile.glint --ewas --wilc --pheno y1

will run EWAS using the Wilcoxon rank-sum test.



|
|

Linear mixed model (LMM)
^^^^^^^^^^^^^^^^^^^^^^^^

.. _--lmm:

**--lmm**

Performs EWAS on the data using linear mixed model (LMM). This is an implementation of the FaST-LMM algorithm by Lippert et al. [1]_

The output file named *results.glint.lmm.txt** includes a list of the sites, sorted by their association p-value. The output file includes the following columns:  ID (CpG identifiers), chromosome (chromosome number of the site), MAPINFO (position of the site in the genome), p-value, q-value, intercept , V1 (coefficient of the first covariate),..., Vn (coefficient of the last covaraite, beta (the coefficient of the site under test), statistic (the test statistic), sigma-e (an estimate of sigma_e), sigma-g (an estimate of sigma_g), UCSC_RefGene_Name (name of the gene that is closest to this site), Relation_to_UCSC_CpG_Island (category)


**--kinship**

The kinship matrix for modelling the inter-individual similarity in the data that is required for the LMM. glint allows two options:

- User-supplied kinship - users can suplly a text file with samples by samples kinship matrix (with no row or column headers). 
- *refactor* - the ReFACTor algorithm can be used for constructing the kinship matrix. If this option is used then ReFACTor is executed for selecting the top informative sites in the data. The kinship matrix is then constructed by calculatign the empirical covariance matrix of the samples based on the selected sites.

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno y1 --kinship kinship.txt

will run EWAS using LMM with the kinship matrix specified in the *kinship.txt* file. Alternatively:

	glint.py --datafile datafile.glint --ewas --lmm --pheno y1 --kinship refactor --k 6

will use the ReFACTor algorithm for constructing the kinship matrix (where 6 is the number of assumed cell types, see the argument `--k`_ for more details).


.. note:: If the *refactor* option is used then all of the arguments available with the `--refactor`_ argument are also available here.



**--reml**

Allows to indicate whether rstricted maximum likelihood estimation (REML) or maximum likelihood estimation (ML) should be used. The default value is 1 (REML). Alternatively, ML can be selected usign the value 0.

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno y1 --kinship kinship.txt --reml 0

will perform EWAS on the data using LMM with ML estimation.


**--norm**

This argument normalizes the covariates (if supplied) before fitting the LMM.

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno y1 --covar c1 c2 c3 --norm

will perform EWAS on the data using LMM after normalizing the covariates c1, c2 and c3.



**--oneld**

This argument allows to fit the log delta parameter in the Fast-LMM model only once (instead for each site separately).

For example::

	glint.py --datafile datafile.glint --ewas --lmm --pheno y1 --oneld

will perform EWAS on the data using LMM with a single value of log detla.



.. _--phenofile: input.html#phenofile

.. _--covarfile: input.html#covarfile

.. _--plot: plots.html#plot

.. _--k: tissueheterogeneity.html#k

.. _--refactor: tissueheterogeneity.html#refactor

.. _glint files: input.html#glint-files


.. [1] Lippert, Christoph, Jennifer Listgarten, Ying Liu, Carl M. Kadie, Robert I. Davidson, and David Heckerman. "FaST linear mixed models for genome-wide association studies." Nature methods 8, no. 10 (2011): 833-835.
