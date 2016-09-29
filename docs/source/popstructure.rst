
Inferring population structure
==============================

Here we provide an implementation of the Epistructure algorithm by Rahmani et al. [1]_ for inferring population structure from methylation without the need for genotyping data. The Epistructure algorithm calculates components that are correlated with the ancestry information of the samples in the data and therefore can be added as covariates in a downstread analysis.


Epistructure
^^^^^^^^^^^^

**--epi:**

Computes the Epistructure components and generates a file named epistructure.pcs.txt with the output.

.. note:: Add "--out filename" to change the default output name.

.. note:: Use --epi together with --gsave in order to generate a new glint file with the computed Epistructure components.

TODO: mention that people should not remove the polymorphic CpGs prior to running Epistructure

**--covar:**

Selects covariates to use in the calculating of the Epistructure components, according to the covaraite names (the headers provided in the covariates file using --covarfile). Considering dominant batch effects and strong genome-wide effectors such as cell type composition is expected to improve the correlation of the Epistructure components with the cell type composition.

For example::

	glint.py --datafile datafile.glint --epi --covar c1,c2,c3

will compute the Epistructure components while using covaraites c1, c2 and c3.


.. note:: Use the argument --covarfile in order to provide covaraites.



--savepcs

Selectes the number of Epistructure components to output (default is 1).

For example::

	glint.py --datafile datafile.glint --epi --savepcs 2

will compute and output the first two Epistructure components.





.. [1] Rahmani, Elior, Liat Shenhav, Regev Schweiger, Paul Yousefi, Karen Huen, Brenda Eskenazi, Celeste Eng et al. "Genome-wide methylation data mirror ancestry information." bioRxiv (2016): 066340.

