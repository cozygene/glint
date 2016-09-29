
Tissue heterogeneity
====================

Tissue heterogenity is known to be a major source of variation and therefore a potential confounder in EWAS. Here we provide an implementation of the ReFACTor algorithm by Rahmani et al. [1]_ for inferring cell compostion.
The ReFACTor algorithm calculates components that are correlated with the cell type composition of the samples in the data and therefore can be added as covariates in a downstread analysis. Note that ReFACTor calculates linear transformations of the cell type composition rather than estiamtes of absolute cell proportion values. 


ReFACTor
^^^^^^^^

**--refactor:**

Computes the ReFACTor components and generates two files:

- refactor.components.txt - a text file with the ReFACTor components.
- refactor.rankedlist.txt - a text file with the CpGs in the data matrix ranked by their level of information according to ReFACTor.

.. note:: Add "--out filename" to change the default output name.

.. note:: Use --refactor together with --gsave in order to generate a new glint file with the computed ReFACTor components.



**--k:**

The assumed number of cell types. This is the only required argument for ReFACTor.

For example::

	glint.py --datafile datafile.txt --refactor --k 6

will compute the ReFACTor components under the assumption of 6 cell types in the data.



**--covar:**

Selects covariates to use in the feature selection step of ReFACTor, according to the covaraite names (the headers provided in the covariates file using --covarfile). Considering genome-wide effectors such as batch information, gender, age, ancestry etc. is expected improve the correlation of the ReFACTor components with the cell type composition.

For example::

	glint.py --datafile datafile.glint --refactor --k 6 --covar c1,c2,c3

will compute the ReFACTor components while using covaraites c1, c2 and c3.


.. note:: Use the argument --covarfile in order to provide covaraites.



**--stdth:**

Excludes sites with standard deviation lower than a specified value. The default value is 0.02.

For example::

	glint.py --datafile datafile.txt --refactor --k 6 --stdth 0.01

will remove all sites with standard deviation lower than 0.01 before computing the ReFACTor components 



**--t:**

The number of sites to use for computing the ReFACTor components. The default value is 500.

For example::

	glint.py --datafile datafile.txt --refactor --k 6 --t 1000

will compute the ReFACTor components using the top 1000 most informative sites according to ReFACTor.



**--numcomp:**

The number of ReFACTor components to output. The default value is the same argument given with --k.

For example::

	glint.py --datafile datafile.txt --refactor --k 6 --numcomp 10

will output the first 10 ReFACTor components.

.. note:: While the argument --k is the assumed number of cell types in the data, --numcomp allows to output more ReFACTor componetns than the number specified by --k (in some cases the number of assumed cell types k may be captured by more than k ReFACTor components).



**--fs:**

The type of feature selection procedure to perform in the feature selection step of ReFACTor.

- normal (default) - the standard feature selection as described in the ReFACTor paper.
- controls - the standard ReFACTor feature selection but based on the control samples only. This option requires the phenotype to be binary (case / control; the controls are assume to be coded as '0'). This option is especially preferable in case where many sites are expected to be assocaited with the phenotype of interest.
- phenotype - a continuous version of the "controls" feature selecction. This feature selection uses the standard ReFACTor feature selection after adjusting the data for the phenotype of interest, in attempt to avoid capturing true signal of the phenotype that is orthogonal to the cell type composition information in the data. This option is especially preferable in case where many sites are expected to be assocaited with the phenotype of interest.

For example::

	glint.py --datafile datafile.txt --refactor --k 6 --fs controls

will compute the ReFACTor components using the "controls" feature selection.


.. note:: If multiple phenotypes were added to the glint file then the --pheno argument should be provided with the name of the phenotypes of interest.



.. [1] Rahmani, Elior, Noah Zaitlen, Yael Baran, Celeste Eng, Donglei Hu, Joshua Galanter, Sam Oh et al. "Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies." Nature methods 13, no. 5 (2016): 443-445.





