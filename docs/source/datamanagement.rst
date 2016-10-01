
Data management
===============

The following section describes arguments that allow to perform basic data management and quality control procedures on data files.

.. note:: Data management commands applied to data do not change the input files. In order to save the changes use the `--gsave`_ or `--save`_ commands.

.. _--gsave: input.html#glint-files

.. _--save: input.html#glint-files


Outliers detection
^^^^^^^^^^^^^^^^^^

**--maxpcstd**


Filters outlier samples using PCA. This argument removes samples with principal components (PCs) above or below a specified number of standard deviations.

For example::

	glint.py --datafile datafile.glint --maxpcstd 1 3

will remove all samples having extreme values of more than 3 standard deviations in their PC 1.

.. note:: Use the --plotpcs argument described under "plots" for plotting the first several PCs of the samples in order to determine whether outliers exist in your data.

.. note:: You can remove outliers based on more than one PC at the same time. For example, for indicating 3 SDs as the maximum level allowed for both PC 1 and PC2 use: *--maxpcstd 1 3 --maxpcstd 2 3*.

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


**--rmxy**


**--rmns**


**--rmpoly**



