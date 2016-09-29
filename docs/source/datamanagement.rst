
Data management
===============

Note that the data management commands are applied to tasks that are added to the command.. the data management commands are then not saved unless --gsave is used

glint files
^^^^^^^^^^^

**--gsave**

Saves a glint binary format of the methylation data. If a phenotype file (--pheno) and covaraites file (--covar) are provided then their information is also contained in the glint file. In addition to the .glint file generated, this argument also generates two additional files:

- datafile.samples.txt - contains the sample identifiers of the samples in the data and the phenotypes and covaraites for each sample.
- datafile.sites.txt - contains the CpGs identifiers of the sites in the data and additional information for each CpG.

For example::

	glint.py --datafile datafile.txt --gsave tutorial_datafile

will save a glint file named tutorial_datafile.glint.



**--save**

Saves the methylation data matrix into a text file.

For example::

	glint.py --datafile datafile.glint --save filename

will save the methylation matrix in filename into a new text file named filename.


Outliers detection
^^^^^^^^^^^^^^^^^^

**--maxpcstd**


Filters outlier samples using PCA. This argument removes samples with principal components (PCs) above or below a specified number of standard deviations.

For example::

	glint.py --datafile datafile.glint --maxpcstd 1 3

will remove all samples having extreme values of more than 3 standard deviations in their PC 1.

.. note:: Use the --plotpcs argument described under "plots" for plotting the first several PCs of the samples in order to determine whether outliers exist in your data.



Data filtering
^^^^^^^^^^^^^^

**--include**

Allows to filter sites. This argument gets a text file with a list of CpG identifiers and considers only these sites in the subsequent analysis.

For example::

	glint.py --datafile datafile.glint --include list.txt

will exclude all the sites not indicated in the list.txt file.



**--exclude**

Allows to filter sites. This argument gets a text file with a list of CpG identifiers and considers all sites except for those in the list in the subsequent analysis.

For example::

	glint.py --datafile datafile.glint --exclude list.txt

will exclude all the sites indicated in the list.txt file.



**--keep**

Allows to filter sampels. This argument gets a text file with a list of sample identifiers and considers only these samples in the subsequent analysis.

For example::

	glint.py --datafile datafile.glint --keep list.txt

will remove all the samples not indicated in the list.txt file.



**--remove**

Allows to filter sampels. This argument gets a text file with a list of sample identifiers and considers all samples except for those in the list in the subsequent analysis.

For example::

	glint.py --datafile datafile.glint --exclude list.txt

will remove all the samples indicated in the list.txt file.



**--minmean**

Filters sites by their mean methylation levels. This argument removes all sites with mean methylation level below the specified value.

For example::

	glint.py --datafile datafile.glint --minmean 0.2

will remove all the sites with mean methylation level bellow 0.2.




**--maxmean**

Filters sites by their mean methylation levels. This argument removes all sites with mean methylation level above the specified value.

For example::

	glint.py --datafile datafile.glint --maxmean 0.8

will remove all the sites with mean methylation level above 0.8.


