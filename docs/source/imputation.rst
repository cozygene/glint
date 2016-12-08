

Imputation from genotype data
=============================

GLINT allows to impute methylation levels from genotype data. Assuming genotype data are available for your samples, GLINT can impute the methylation levels of some methylation sites.

The imputation is based on the EPISTRUCTURE paper by Rahmani et al. [1]_. As described in the paper, a linear model was fitted for each methylation site from cis-located SNPs in a large dataset for which both methylation and genotype data were available. A score was then defined for each site, based on the squared linear correlation of the model. Here, we use the coefficients of these linear models in order to predict methylation levels based on the SNPs. Sites having higher scores are expected to be predicted more accurately compared with sites having lower scores.

.. note:: The example commands described bellow assume that the user generated `GLINT files`_ with covariates file and phenotypes file.

.. note:: We suggest users to apply a full GWAS quality control pipeline on their genotype data before imputing metylation levels, e.g. using PLINK_.

.. note:: The linear models for imputation were fitted based on European individuals.

.. note:: A/T and C/G SNPs are not used for the imputation, in order to avoid strand misspecification due to different genotyping platforms.

.. _--impute:

**--impute:**

Imputes methylation levels from genotype data. This argument requires using 3 additional arguments for loading the genotype data in an Eigenstrat format.
You can find here_ details about the EIGENSTRAT format.


For example::

	python glint.py --impute --snp genotypes.snp --geno genotypes.geno --ind genotypes.ind

will generate GLINT files with imputed methylation data for a group of methylation sites, based on the *.snp*, *.geno* and *.ind* EIGENSTRAT files of the genotype data.

.. note:: --impute will automatically generate GLINT files with the imputed methylation levels, therefore there is no need to add the --gsave argument.

.. note:: Polymorphic CpGs according to Chen et al. [2]_ are not imputed. In addition, methylation levels are imputed only for CpGs for which at least one of their predicting SNPs exists in the provided genotype data.

.. note:: Use `--out`_ in order to change the default output name.


.. _--score:

**--score**

Controls the number of methylation sites to impute. This argument specifies the minimal score required for a site in order to get imputed. The deafult value is 0.5.


For example::

	python glint.py --impute --snp genotypes.snp --geno genotypes.geno --ind genotypes.ind --score 0.4

will impute methylation levels for every site with a score greater than 0.4.


.. _--maxmiss”:

**--maxmiss**



Controls the fraction of missing values allowed for SNPs in the EIGENSTRAT input files. SNPs with fraction of missing values above the specified value will be ignored and will not be used in the imputation. The deafult value is 0.03.

For example::

	python glint.py --impute --snp genotypes.snp --geno genotypes.geno --ind genotypes.ind --maxmiss 0.01

will impute methylation levels for the samples in the data using SNPs with no more than 1% of missing values.




.. _--out: input.html#out

.. _GLINT files: input.html#glint-files

.. _PLINK: http://pngu.mgh.harvard.edu/~purcell/plink/


.. _here: http://genepath.med.harvard.edu/~reich/InputFileFormats.htm

.. [1] Rahmani, Elior, Liat Shenhav, Regev Schweiger, Paul Yousefi, Karen Huen, Brenda Eskenazi, Celeste Eng et al. "Genome-wide methylation data mirror ancestry information." bioRxiv (2016): 066340.

.. [2] Chen, Yi-an, Mathieu Lemire, Sanaa Choufani, Darci T. Butcher, Daria Grafodatskaya, Brent W. Zanke, Steven Gallinger, Thomas J. Hudson, and Rosanna Weksberg. "Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray." Epigenetics 8, no. 2 (2013): 203-209.

