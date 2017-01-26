


Tissue heterogeneity
====================


When methylation data are collected from an heterogeneous source (e.g. whole blood) the cell type compositions of samples in the data are known to be a major source of variation, and therefore a potential confounder in EWAS. Here we provide an implementation of the ReFACTor algorithm by Rahmani et al. for inferring cell type compostion information in methylation data [1]_. In addition, we provide an implementation of the reference-based method by Houseman et al. [2]_, allowing to estimate cell counts (cell proportion) of the samples in the data.

.. note:: The example commands described bellow assume that the user generated `GLINT files`_ with covariates file and phenotypes file.


|
|

ReFACTor
^^^^^^^^

The ReFACTor algorithm calculates components that are correlated with the cell type composition of the samples in the data by applying an unsupervised feature selection step followed by principal component analysis (PCA). These components can be added as covariates in a downstread analysis in order to account for the tissue heterogeneity. Note that ReFACTor calculates linear transformations of the cell type composition rather than estimates of absolute cell count values.

.. _--refactor:

**--refactor**

Computes the ReFACTor components and generates two output files:

- *refactor.components.txt* - a text file with the ReFACTor components.
- *refactor.rankedlist.txt* - a text file with the CpGs in the data matrix ranked by their level of information according to ReFACTor.


.. note:: For best performance of ReFACTor we recommend adding known gneome-wide effectors as covariates using the `--covar`_ argument.

.. note:: GLINT computes the ReFACTor components while automatically ignoring X and Y chromosomes sites, polymorphic sites and cross-reactive sites according to Chen et al. [3]_ for 450K array data and according to McCartney et al. [4]_ for 850K (EPIC) array data.

.. note:: Use `--out`_ in order to change the default output name.

.. note:: Use `--refactor`_ together with `--gsave`_ in order to generate a new version of GLINT files with the computed ReFACTor components (these will be included in the *datafile.samples.txt* file).




.. _--k:

**--k**

The assumed number of cell types in the data. This is the only required argument for ReFACTor.

For example::

	python glint.py --datafile datafile.glint --refactor --k 6

will compute the ReFACTor components under the assumption of 6 cell types in the data.


.. _--covar:

**--covar**

Selects covariates to use in the feature selection step of ReFACTor. Considering genome-wide effectors such as batch information, gender, age, ancestry etc. is expected to improve the correlation of the ReFACTor components with the cell type composition.

For example::

	python glint.py --datafile datafile.glint --refactor --k 6 --covar c1 c2 c3

will compute the ReFACTor components while accounting for the covariates c1, c2 and c3. The names of the covariates are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `GLINT files`_.

.. note:: Use the argument `--covarfile`_ in order to provide covariates that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.


.. _--stdth:

**--stdth**

Excludes sites with standard deviation lower than a specified value. This argument results in a speed-up of ReFACTor while ignoring sites with low variability that are less likely to contain cell composition information. The default value is 0.02.

For example::

	python glint.py --datafile datafile.glint --refactor --k 6 --stdth 0.01

will remove all sites with standard deviation lower than 0.01 before computing the ReFACTor components.


.. _--t:

**--t**

The number of sites to use for computing the ReFACTor components. The default value is 500.

For example::

	python glint.py --datafile datafile.glint --refactor --k 6 --t 1000

will compute the ReFACTor components using the top 1000 most informative sites according to ReFACTor.


.. _--numcomp:

**--numcomp**

The number of ReFACTor components to output. The default value is the same value given with the `--k`_ argument.


For example::

	python glint.py --datafile datafile.glint --refactor --k 6 --numcomp 10

will output the first 10 ReFACTor components.

.. note:: While the argument `--k`_ is the assumed number of cell types in the data, `--numcomp`_ allows to output more or less ReFACTor components than the number specified by `--k`_ (in some cases the number of assumed cell types k may be captured by more than k ReFACTor components).


.. _--fs:

**--fs**

The type of feature selection procedure to perform in the feature selection step of ReFACTor.

- *normal* (default) - the standard feature selection as described in the ReFACTor paper.
- *controls* - the standard ReFACTor feature selection but based on the control samples only. This option requires the phenotype to be binary (case / control; the controls are assume to be coded as '0'). This option is especially favourable in case where many sites are expected to be assocaited with the phenotype of interest.
- *phenotype* - a continuous version of the *controls* feature selection. This feature selection uses the standard ReFACTor feature selection after adjusting the data for the phenotype of interest, in attempt to avoid capturing true signal of the phenotype that is independent in the cell type composition information in the data. This option is especially favourable in case where many sites are expected to be assocaited with the phenotype of interest.

For example::

	python glint.py --datafile datafile.glint --refactor --k 6 --fs controls --pheno y1

will compute the ReFACTor components using the *controls* feature selection based on the phenotype y1. The names of the phenotypes are defined by the headers in the *datafile.samples.txt* file associated with the *datafile.glint*. For more details see `GLINT files`_.

.. note:: Use the argument `--phenofile`_ in order to provide phenotypes that were not included in the *datafile.glint* file or in case where a textual version of the data is used rather than a *.glint* file.



Houseman
^^^^^^^^

The algorithm by Houseman et al. is a reference-based method for calculating cell count estimates. This method requires reference data of cell type specific mean methylation levels of sorted cell types from the studied tissue. The default reference data is based on whole-blood data by Reinius et al. [5]_, according to the eature selection proposed by Koestler et al [6]_.

.. note:: Reference data currently exist for 7 leukocyte cell types only.



.. _--houseman:

**--houseman**

Computes cell count estimates according to the algorithm by Houseman et al. and generates an output file titled *houseman_estimates.txt*, containing cell count estimates.

For example::

	python glint.py --datafile datafile.glint --houseman

will compute cell count estimates.


.. note:: Use `--out`_ in order to change the default output name.

.. note:: Use `--houseman`_ together with `--gsave`_ in order to generate a new version of GLINT files with the computed cell count estimates (these will be included in the *datafile.samples.txt* file).


.. _--reference:

**--reference**

Allows to include user-supplied reference data. This argument gets path to a file containing sites by cell types matrix of mean methylation levels (for each methylation site in each cell type). The first row should include "ID" followed by cell type names and the first column should include CpG identifiers. The default reference in GLINT contains 7 leukocyte cell types. The file can be either tab-delimited, comma-delimited or space-delimited.


For example::

	python glint.py --datafile datafile.glint --houseman --reference reference.txt

will compute cell count estimates using the reference data in *reference.txt*.








.. _--gsave: input.html#gsave

.. _--out: input.html#out

.. _--covarfile: input.html#covarfile

.. _--phenofile: input.html#phenofile

.. _GLINT files: input.html#glint-files





.. [1] Rahmani, Elior, Noah Zaitlen, Yael Baran, Celeste Eng, Donglei Hu, Joshua Galanter, Sam Oh et al. "Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies." Nature methods 13, no. 5 (2016): 443-445.

.. [2] Houseman, Eugene Andres, William P. Accomando, Devin C. Koestler, Brock C. Christensen, Carmen J. Marsit, Heather H. Nelson, John K. Wiencke, and Karl T. Kelsey. "DNA methylation arrays as surrogate measures of cell mixture distribution." BMC bioinformatics 13, no. 1 (2012): 1.

.. [3] Chen, Yi-an, Mathieu Lemire, Sanaa Choufani, Darci T. Butcher, Daria Grafodatskaya, Brent W. Zanke, Steven Gallinger, Thomas J. Hudson, and Rosanna Weksberg. "Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray." Epigenetics 8, no. 2 (2013): 203-209.

.. [4] McCartney, Daniel L., Rosie M. Walker, Stewart W. Morris, Andrew M. McIntosh, David J. Porteous, and Kathryn L. Evans. "Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip." Genomics Data 9 (2016): 22-24.

.. [5] Reinius, Lovisa E., Nathalie Acevedo, Maaike Joerink, Göran Pershagen, Sven-Erik Dahlén, Dario Greco, Cilla Söderhäll, Annika Scheynius, and Juha Kere. "Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility." PloS one 7, no. 7 (2012): e41361.

.. [6] Koestler, Devin C., Meaghan J. Jones, Joseph Usset, Brock C. Christensen, Rondi A. Butler, Michael S. Kobor, John K. Wiencke, and Karl T. Kelsey. "Improving cell mixture deconvolution by id entifying o ptimal DNA methylation l ibraries (IDOL)." BMC bioinformatics 17, no. 1 (2016): 1.
