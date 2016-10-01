# glint
An easy and efficient command-line tool for the analysis of DNA methylation data and conducting EWAS.  
Our software does not require any programming proficiency.

Some of the abaliable features:
 - Data managment
 - Correction of cell type composition (ReFACTor)
 - Ancestry estimation
 - A novel method for methylation imputation 
 - EWAS tests: linear regressin, logistic regression, LMM, Wilcoxon
 - Generation of publication-quality plots

For more details see "Documentation"

### Download and Installation

1. Download the latest release from <a href="put the link here todo" target="_blank"> todo!edit this!here</a>.
2. Install <a href="https://www.continuum.io/downloads" target="_blank">Anaconda Python version 2.7</a> (automatically includes all required dependencies).

  If you already have Python installed and do not want to install Anaconda Python, skip to step 3.

3. Install todoPUTNAMEHERE by running the provided install.py file, run: ```python install.py```

For more details see "Dependencies".
  
### Quick example
This command will conduct EWAS on a metylation data and a phenotype (files 'data.txt' and 'phenotype.txt').
Fisrt it will remove sites low variance sites, then it will run association test (linear regression) for each site on a phenotype. Finally it will plot the results and will generate a publication-quality qq-plot and a Manhattan plot.
```
python glinttodo.py --datafile data.txt --minstd 0.02 --ewas --linreg --phenofile phenotype.txt --pheno age --plot --qqplot --manhhattan
```
*In order to quick learn more options, try this <a href="todo add link to tutorial" target="_blank">tutorial</a>*

### Documentation
Detailed documentation explaining all the features can be found <a href="todo add link to docs" target="_blank">here</a>.  
We also supply a quick  <a href="todo add link to tutorial" target="_blank">tutorial</a> that will walk you through the basic options avaliable

### Dependencies

This release of todoPUTNAMEHERE was implemented for Python 2.7 and has the following dependencies:

    numpy
    scipy
    sklearn
    pandas
    matplotlib
    statsmodels
    

We recommend installing <a href="https://www.continuum.io/downloads" target="_blank">Anaconda Python version 2.7</a>, which already includes all necessary dependencies.

If you already have Python installed and do not want to install Anaconda Python, run "install.py" script (found in the "python" folder):
```
python install.py
```
The script automatically installs missing dependencies that are required for todoPUTNAMEHERE. Note that in some environments the script may fail to install some of the dependencies, in which case you will need to manually install them.

### Citing todoPUTNAMEHERE
elior..todo


### Authors

This software was developed by Reut Yedidim, Omer Weissbrod  and Elior Rahmani.

For any question and for reporting bugs please send an email to Elior Rahmani at: elior.rahmani@gmail.com
