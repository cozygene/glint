
Imputation from genotype data
=============================

glint allows to impute methylation levels from genotyping data. Assuming genotype data are available for your samples, glint can impute the methylation levels in some methylation sites using a linear model. The number of imputed methylation sites depends on the .. and on the quality of imputation as was.. (see the --score argument for more details).

The imputation is based on the epistructure paper - cite... [2]_.


**--impute:**

Imputes methylation levels from genotype data. This argument requires using 3 additional arguments --snp, --in and --geno. As described bellow, these arguments are used to load the genotype data files in an EIGENSTRAT format (for details about the EIGENSTRAT format see [1]_).

For example::

	glint.py --impute --snp genotypes.snp --geno genotypes.geno --ind genotypes.ind

will generate glint files with imputed methylation data for a group of methylation sites.

.. note:: --impute will automatically generate glint files with the imputed methylation levels, therefore there is no need to add the --save argument.

.. note:: Add "--out filename" to change the default output name.


**score**


**maxmiss**






.. [1] http://genepath.med.harvard.edu/~reich/InputFileFormats.htm

.. [2] Rahmani, Elior, Liat Shenhav, Regev Schweiger, Paul Yousefi, Karen Huen, Brenda Eskenazi, Celeste Eng et al. "Genome-wide methylation data mirror ancestry information." bioRxiv (2016): 066340.
