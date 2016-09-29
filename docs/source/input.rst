
Input files
===========

The following arguments are used to include data files. Example files can be downloaded here_. 

.. _here: blank


**--datafile:**	

Path to a tab-delimited file of beta-normalized sites by samples methylation matrix. The first row should include "ID" followed by sample identifiers and the first column should include CpG identifiers. Note that a data file must be provided with every glint command. For example file, see tutorial_datafile.txt.

For example::

	glint.py --datafile tutorial_datafile.txt

will load the methylation data matrix in the tutorial_datafile.txt file.

.. note:: For users having Rdata files with methylation data we provide a script for generating files in the required format - see "Convert RData file into glint format" at the end of this page.

.. note:: In order to speed up glint, we recommend working with a binary version of the data files. See **--gsave** under "Data management" for more details.



**--covarfile**

Path to a tab-delimited file of samples by covariates matrix. The first row can include "ID" followed by headers, one for each covariate, and the first column should include sample identifiers. If a row of headers is not provided then glint will automatically generate a name for each covariate. For example file, see tutorial_covariates.txt.

For example::

	glint.py --datafile tutorial_datafile.txt --covarfile tutorial_covariates.txt

will load the covariate data matrix in the tutorial_covariates.txt file.

TODO: more than one file can be given..

**--phenofile**

Path to a tab-delimited file of samples by phenotypes matrix. The first row can include "ID" followed by headers, one for each phenotype, and the first column should include sample identifiers. If a row of headers is not provided then glint will automatically generate a name for each phenotype. For example file, see tutorial_phenotypes.txt.

For example::

	glint.py --datafile tutorial_datafile.txt --phenofile tutorial_phenotypes.txt

will load the phenotypes data matrix in the tutorial_phenotypes.txt file.

TODO: more than one file can be given..


Convert RData file into glint format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**convertToGlintInput.R:**

We provide this R script for users having methylation data matrix in RData format. This script gets as an input RData file with sites by samples methylation data matrix saved as a data frame or matrix variable. The data matrix should incude CpGs identifiers as volumn names and sample identifiers as row names. The script takes as an input the RData file name, and optionally can take two additional arguments:

- varname - if more than a single data frame / matrix variable exists in the RData file then the name of the methylation data matrix variable should be provided. If this argument is not provided then the script automatically attemps to find data frame or a matrix variable.
- transpose - if the methylation data matrix is formatted as samples by sites rather than sites by samples then providing this argument with the value 'true' will transpose the data matrix.

For example::

	Rscript convertToGlintInput.R datafile.RData X

will save a tab-delimited text file containing sites by samples methylation data matrix as appear in the variable X in datafile.RData. This file can be then provided as an input to gilnt (using --datafile).
Alternatively::

	Rscript convertToGlintInput.R datafile.RData X true

will assume that the information in the variable X is formatted as samples by sites and therefore should be transposed.


