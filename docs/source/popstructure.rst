

Inferring population structure
==============================

Here we provide an implementation of the EPISTRUCTURE algorithm by Rahmani et al. [1]_ for inferring population structure from methylation without the need for genotyping data. The EPISTRUCTURE algorithm calculates components that are correlated with the ancestry information of the samples in the data by applying principal component analysis (PCA) on a set of pre-defined list of sites that were found to capture a high level of genetic information. The EPISTRUCTURE components can be added as covariates in a downstream analysis.

.. note:: The example commands described bellow assume that the user generated `GLINT files`_ with covariates file and phenotypes file.

.. note:: The reference list of sites were found based on European individuals.


|
|

EPISTRUCTURE
^^^^^^^^^^^^

.. _--epi:

**--epi:**

Computes the EPISTRUCTURE components and generates a file titled *epistructure.pcs.txt* with the output.

For example::

	python glint.py --datafile datafile.glint --epi

will compute the EPISTRUCTURE components of the data.


.. note:: EPISTRUCTURE leverages polymorphic sites in order to capture the genetic and therefore the ancesty information in the data better. Therefore, we recommend to avoid removing polymorphic sites (`--rmpoly`_) before applying EPISTRUCTURE.

.. note:: Use `--epi`_ together with `--gsave`_ in order to generate a new version of GLINT files with the computed EPISTRUCTURE components (these will be included in the *datafile.samples.txt* file).

.. note:: Use `--out`_ in order to change the default output name.


.. _--covar:

**--covar:**

Selects covariates to use in the calculation of the EPISTRUCTURE components. Considering highly dominant genome-wide effectors such as cell type composition (in case of heterogeneous tissue) is expected to improve the correlation of the EPISTRUCTURE components with the cell type composition.

For example::

	python glint.py --datafile datafile.glint --epi --covar c1 c2 c3

will compute the EPISTRUCTURE components while accounting for the covariates c1, c2 and c3. The names of the covariates are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `GLINT files`_.

.. note:: Use the argument `--covarfile`_ in order to provide covariates that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.




.. _--savepcs:

--savepcs

Selectes the number of EPISTRUCTURE components to output (default is 1).

For example::

	python glint.py --datafile datafile.glint --epi --savepcs 2

will compute the first two EPISTRUCTURE components of the data.



.. _--covarfile: input.html#covarfile

.. _--gsave: input.html#gsave

.. _--out: out.html#out

.. _--rmpoly: datamanagement.html#rmpoly

.. _GLINT files: input.html#glint-files



.. [1] Rahmani, Elior, Liat Shenhav, Regev Schweiger, Paul Yousefi, Karen Huen, Brenda Eskenazi, Celeste Eng et al. "Genome-wide methylation data mirror ancestry information." bioRxiv (2016): 066340.

