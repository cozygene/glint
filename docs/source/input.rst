
Input
=====

The following section describes arguments that allow to provide data files for glint, and the *--gsave* argument that allows to save and work with a binary version of the data (glint files) for gaining computation speed-up.
Note that glint does not work on raw data, but rather assumes data are given after raw data preprocessing. 



Input files
^^^^^^^^^^^

The `tutorial files`_ can be used as example files.

.. _tutorial files: blank

**--datafile:**	

Path to a tab-delimited file of sites by samples matrix of methylation levels. The first row should include "ID" followed by sample identifiers and the first column should include CpG identifiers. 

For example::

	glint.py --datafile tutorial_datafile.txt

will load the methylation data matrix in the *datafile.txt* file. See the tutorial *datafile.txt* file as an example file.

.. note:: For users having *.Rdata* files with methylation data we provide a script for generating files in the required format - see `Convert RData file into glint format`_ for more details.

.. note:: In order to speed-up glint, we recommend working with a binary version of the data. See `--gsave`_ for more details.



**--covarfile**

Path to a tab-delimited file of samples by covariates matrix. The first row may be a row of headers - "ID" followed by the names of the covariates, and the first column should include sample identifiers. If a row of headers is not provided then glint will automatically generate a name for each covariate.

For example::

	glint.py --datafile datafile.txt --covarfile covariates.txt

will provide the covariates matrix in the *covariates.txt* file. See the tutorial *covariates.txt* file as an example file.

.. note:: More than one covariates file can be provided, e.g. *--covarfile covariates1.txt covariates2.txt*.


**--phenofile**

Path to a tab-delimited file of samples by phenotypes matrix. The first row may be a row of headers - "ID" followed by the names of the phenotypes, and the first column should include sample identifiers. If a row of headers is not provided then glint will automatically generate a name for each phenotype.

For example::

	glint.py --datafile datafile.txt --phenofile phenotypes.txt

will provide the phenotypes matrix in the *phenotypes.txt* file. See the tutorial *phenotypes.txt* file as an example file.


.. note:: More than one phenotypes file can be provided, e.g. *--phenofile phenotypes1.txt phenotypes2.txt*.


|
|

.. _--gsave:

glint files
^^^^^^^^^^^

**--gsave**

Saves glint files, including a binary version of the methylation data (*.glint* file) and two additional files:

- *datafile.sites.txt* - contains the CpG identifiers of the sites in the data and additional information for each CpG: chromosome, position, nearest gene and genomic category.

- *datafile.samples.txt* - contains the sample identifiers of the samples in the data. If *--covarfile* and *--phenofile* are used then this file also includes the phenotypes and covaraites for each sample.

For example::

	glint.py --datafile datafile.txt --gsave datafile

will save a binary data file titled *datafile.glint* and two additional files titled *datafile.samples.txt* and *datafile.sites.txt*. The following command:

::

	glint.py --datafile datafile.txt --gsave datafile --covarfile covariates.txt --phenofile phenotypes.txt

will also include the covariates and phenotypes information found in the *covariates.txt* and *phenotypes.txt* files in the *datafile.samples.txt* file.



**--save**

Allows to save a textual version of the data in a binary *.glint* file.

For example::

	glint.py --datafile datafile.glint --save datafile.txt

will create a file titled *datafile.txt* with a textual version of the methylation matrix in *datafile.glint*.

.. note:: *--save* can be also used to save a new version of textual format of previous textual files (i.e. *--save* is not restricted to get *.glint* file as an input).


|
|

.. _Convert RData file into glint format:

Convert R file to glint format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**convertToGlintInput.R:**

We provide this R script for users having methylation data matrix in *.RData* format. This script gets as an input *.RData* file with sites by samples methylation data matrix saved as a data frame or a matrix variable with CpGs identifiers as row names and sample identifiers as column names. In addition to the *.RData* file name, the script optionally can take two additional arguments:

- varname - if more than a single data frame / matrix variable exists in the *.RData* file then the name of the methylation data variable should be provided. If this argument is not provided then the script automatically attemps to find data frame or a matrix variable.
- transpose - if the methylation data matrix is formatted as samples by sites rather than sites by samples then providing this argument with the value 'true' will transpose the data matrix.

For example::

	Rscript convertToGlintInput.R datafile.RData X

will save a tab-delimited text file containing sites by samples methylation data matrix as appear in the variable X that is saved in the *datafile.RData* file. The resulted file can be then provided as an input to glint (using --datafile).

|

Alternatively::

	Rscript convertToGlintInput.R datafile.RData X true

will assume that the information in the variable X is formatted as samples by sites and therefore should be transposed.


