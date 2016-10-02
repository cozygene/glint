


Tissue heterogeneity
====================

When methylation data are collected from an heterogeneous source (e.g. whole blood) the cell type compositions of samples in the data are known to be a major source of variation, and therefore a potential confounder in EWAS. Here we provide an implementation of the ReFACTor algorithm by Rahmani et al. for inferring cell type compostion in methylation data [1]_.


.. note:: The example commands described bellow assume that the user generated `glint files`_ with covariates file and phenotypes file.


|
|

ReFACTor
^^^^^^^^

The ReFACTor algorithm calculates components that are correlated with the cell type composition of the samples in the data by applying an unsupervised feature selection step followed by principal component analysis (PCA). These components can be added as covariates in a downstread analysis in order to account for the tissue heterogeneity. Note that ReFACTor calculates linear transformations of the cell type composition rather than estimates of absolute cell proportion values. 

.. _--refactor:

**--refactor**

Computes the ReFACTor components and generates two output files:

- *refactor.components.txt* - a text file with the ReFACTor components.
- *refactor.rankedlist.txt* - a text file with the CpGs in the data matrix ranked by their level of information according to ReFACTor.


.. note:: For best performance of ReFACTor we recommend adding known gneome-wide effectors as covariates using the `--covar`_ argument.

.. note:: Add *--out filename* in order to change the default output name.

.. note:: Use `--refactor`_ together with `--gsave`_ in order to generate a new version of glint files with the computed ReFACTor components (these will be included in the *datafile.samples.txt* file).



.. _--k:

**--k**

The assumed number of cell types in the data. This is the only required argument for ReFACTor.

For example::

	glint.py --datafile datafile.txt --refactor --k 6

will compute the ReFACTor components under the assumption of 6 cell types in the data.


.. _--covar:

**--covar**

Selects covariates to use in the feature selection step of ReFACTor. Considering genome-wide effectors such as batch information, gender, age, ancestry etc. is expected to improve the correlation of the ReFACTor components with the cell type composition.

For example::

	glint.py --datafile datafile.glint --refactor --k 6 --covar c1 c2 c3

will compute the ReFACTor components while accounting for the covariates c1, c2 and c3. The names of the covariates are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `glint files`_.

.. note:: Use the argument `--covarfile`_ in order to provide covariates that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.


.. _--stdth:

**--stdth**

Excludes sites with standard deviation lower than a specified value. This argument results in a speed-up of ReFACTor while ignoring sites with low variability that are less likely to contain cell composition information. The default value is 0.02.

For example::

	glint.py --datafile datafile.glint --refactor --k 6 --stdth 0.01

will remove all sites with standard deviation lower than 0.01 before computing the ReFACTor components.


.. _--t:

**--t**

The number of sites to use for computing the ReFACTor components. The default value is 500.

For example::

	glint.py --datafile datafile.glint --refactor --k 6 --t 1000

will compute the ReFACTor components using the top 1000 most informative sites according to ReFACTor.


.. _--numcomp:

**--numcomp**

The number of ReFACTor components to output. The default value is the same value given with the `--k`_ argument.


For example::

	glint.py --datafile datafile.glint --refactor --k 6 --numcomp 10

will output the first 10 ReFACTor components.

.. note:: While the argument `--k`_ is the assumed number of cell types in the data, `--numcomp`_ allows to output more or less ReFACTor components than the number specified by `--k`_ (in some cases the number of assumed cell types k may be captured by more than k ReFACTor components).


.. _--fs:

**--fs**

The type of feature selection procedure to perform in the feature selection step of ReFACTor.

- *normal* (default) - the standard feature selection as described in the ReFACTor paper.
- *controls* - the standard ReFACTor feature selection but based on the control samples only. This option requires the phenotype to be binary (case / control; the controls are assume to be coded as '0'). This option is especially favourable in case where many sites are expected to be assocaited with the phenotype of interest.
- *phenotype* - a continuous version of the *controls* feature selection. This feature selection uses the standard ReFACTor feature selection after adjusting the data for the phenotype of interest, in attempt to avoid capturing true signal of the phenotype that is independent in the cell type composition information in the data. This option is especially favourable in case where many sites are expected to be assocaited with the phenotype of interest.

For example::

	glint.py --datafile datafile.glint --refactor --k 6 --fs controls --pheno y1

will compute the ReFACTor components using the *controls* feature selection based on the phenotype y1. The names of the phenotypes are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `glint files`_.

.. note:: Use the argument `--phenofile`_ in order to provide phenotypes that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.



.. _--gsave: input.html#gsave

.. _--covarfile: input.html#covarfile

.. _--phenofile: input.html#phenofile

.. _glint files: input.html#glint-files





.. [1] Rahmani, Elior, Noah Zaitlen, Yael Baran, Celeste Eng, Donglei Hu, Joshua Galanter, Sam Oh et al. "Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies." Nature methods 13, no. 5 (2016): 443-445.

