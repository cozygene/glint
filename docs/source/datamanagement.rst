

Data management
===============

The following section describes arguments that allow to perform basic data management and quality control procedures on data.


.. note:: The example commands described bellow assume that the user generated `glint files`_ with covariates file and phenotypes file.

.. note:: Data management commands applied to data do not change the input files. In order to save the changes use the `--gsave`_ or `--txtsave`_ commands.


|
|

Outliers detection
^^^^^^^^^^^^^^^^^^

.. _--maxpcstd:

**--maxpcstd**


Filters outlier samples using PCA. This argument removes samples with principal components (PCs) above or below a specified number of standard deviations.

For example::

	glint.py --datafile datafile.glint --maxpcstd 1 3

will remove all samples having extreme values of PC 1 (more than 3 standard deviations).

.. note:: Use the `--plotpcs`_ argument for plotting the first several PCs of the samples in order to determine whether outliers exist in your data.


.. note:: You can remove outliers based on more than one PC at the same time. For example, for indicating 3 SDs as the maximum level allowed for both PC 1 and PC2 use: *--maxpcstd 1 3 --maxpcstd 2 3*.

|
|

Data filtering
^^^^^^^^^^^^^^

.. _--include:

**--include**

Allows to filter sites. This argument gets a text file with a list of CpG identifiers (one CpG identifier per row) and considers only these sites.

For example::

	glint.py --datafile datafile.glint --include list.txt

will exclude all the sites not indicated in the *list.txt* file.


.. _--exclude:

**--exclude**

Allows to filter sites. This argument gets a text file with a list of CpG identifiers (one CpG identifier per row) and considers all sites except for those in the list.

For example::

	glint.py --datafile datafile.glint --exclude list.txt

will exclude all the sites indicated in the *list.txt* file.


.. _--keep:

**--keep**

Allows to filter sampels. This argument gets a text file with a list of sample identifiers (one sample identifier per row) and considers only these samples.

For example::

	glint.py --datafile datafile.glint --keep list.txt

will remove all the samples not indicated in the *list.txt* file.


.. _--remove:

**--remove**

Allows to filter sampels. This argument gets a text file with a list of sample identifiers (one sample identifier per row) and considers all samples except for those in the list.

For example::

	glint.py --datafile datafile.glint --remove list.txt

will remove all the samples indicated in the *list.txt* file.


.. _--stdth:

**--stdth**

Filters sites by their variance. This argument removes all sites with standard deviation below the specified value. This command can be used to remove nearly-constant sites that are less likely to result in significant associations in EWAS.

For example::

	glint.py --datafile datafile.glint --stdth 0.01

will remove all sites with standard deviation bellow 0.01.


.. _--minmean:

**--minmean**

Filters sites by their mean methylation levels. This argument removes all sites with mean methylation level below the specified value.

For example::

	glint.py --datafile datafile.glint --minmean 0.2

will remove all sites with mean methylation level bellow 0.2.



.. _--maxmean:

**--maxmean**

Filters sites by their mean methylation levels. This argument removes all sites with mean methylation level above the specified value.

For example::

	glint.py --datafile datafile.glint --maxmean 0.8

will remove all sites with mean methylation level above 0.8.



.. _--rmxy:

**--rmxy**

Filters out non-autosomal sites (sites in chromsomes X and Y). This argument assumes that the data were collected using the Illumina 450K array.

For example::

	glint.py --datafile datafile.glint --rmxy

will remove all non-autosomal sites from the data.


.. _--rmns:

**--rmns**

Filters out cross-reactive (non specific) sites according to Chen et al. [1]_. This argument assumes that the data were collected using the Illumina 450K array.

For example::

	glint.py --datafile datafile.glint --rmns

will remove all non specific sites from the data.


.. _--rmpoly:

**--rmpoly**

Filters out polymorphic sites according to Chen et al. [1]_. This argument assumes that the data were collected using the Illumina 450K array.

For example::

	glint.py --datafile datafile.glint --rmpoly

will remove all polymorphic sites from the data.



.. _--gsave: input.html#gsave

.. _--txtsave: input.html#save

.. _--plotpcs: plots.html#plotpcs

.. _glint files: input.html#glint-files


.. [1] Chen, Yi-an, Mathieu Lemire, Sanaa Choufani, Darci T. Butcher, Daria Grafodatskaya, Brent W. Zanke, Steven Gallinger, Thomas J. Hudson, and Rosanna Weksberg. "Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray." Epigenetics 8, no. 2 (2013): 203-209.
