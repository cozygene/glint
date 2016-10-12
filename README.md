# GLINT
An easy and efficient command-line tool for the analysis of DNA methylation data and conducting EWAS.  
Our software does not require any programming proficiency.

The software provides options for data managment, novel methods for correction of cell type composition (ReFACTor), ancestry estimation and methylation imputation. Statistical tests such as linear regressin, logistic regression, LMM and Wilcoxon. Generation of publication-quality plots and more.

For more details see "Documentation"  

### Download and Installation

1. Download the latest release from <a href="put the link here todo" target="_blank"> todo!edit this!here</a>.
2. Install <a href="https://www.continuum.io/downloads" target="_blank">Anaconda Python version 2.7</a> which includes most of our dependencies.  
    - If you already have Python2.7 and don't want to install Anaconda, please see "Dependencies".
3. Install *cvxopt* using Anaconda:  
    **Windows** run ```todo```
    **Linux** run ```sudo `which conda` install -c anaconda cvxopt```
    **MacOS** run ```todo```
    
For more details see "Dependencies".  
  
### Quick example
This command will conduct EWAS on metylation data and a phenotype. First it'll remove the low variance sites, then it'll run linear regression and finally it'll generate a publication-quality qq-plot and Manhattan plot.
```
python glint.py --datafile data.txt --minstd 0.02 --ewas --linreg --phenofile phenotype.txt --pheno age --plot --qqplot --manhhattan
```
**In order to quick learn more options, try this <a href="todo add link to tutorial" target="_blank">tutorial</a>**  

### Documentation
Detailed documentation explaining all the features can be found <a href="todo add link to docs" target="_blank">here</a>.  
We also supply a quick  <a href="todo add link to tutorial" target="_blank">tutorial</a> that will walk you through the basic options avaliable  

### Dependencies

This release of GLINT was implemented for Python 2.7 and has the following dependencies:

    numpy
    scipy
    sklearn
    pandas
    matplotlib
    statsmodels
    cvxopt (not included in Anaconda)
    

We recommend installing <a href="https://www.continuum.io/downloads" target="_blank">Anaconda Python version 2.7</a>, which already includes most of necessary dependencies.  
To install dependency *cvxopt* with Anaconda run:
on **Windows** : ```todo```
on **Linux**: ```sudo `which conda` install -c anaconda cvxopt```
on **MacOS**: ```todo```

If you already have Python installed and do not want to install Anaconda Python, run "install.py" script (found in the "python" folder):
```
python install.py
```
The script automatically installs missing dependencies that are required for GLINT. Note that in some environments the script may fail to install some of the dependencies, in which case you will need to manually install them or  <a href="https://www.continuum.io/downloads" target="_blank">Anaconda</a>.

Note that GLINT firt checkes that all it's dependencies are installed before executing. If it finds a problem it will instruct you.

### Troubleshooting
1. Make sure you run GLINT with Anaconda Python command line, to find it:
  **Windows**
  **Linux** run on command line: ```which conda```.  
    If that command returns nothing than you dont have Anaconda installed, refer to Download and Installation.
    Otherwise, if for example the output of the command is /home/user/anaconda2/bin/conda than Anaconda Python command line is at /home/user/anaconda2/bin/python.
  **MacOs**
2. Make sure you have all dependencies installed.If you installed Anaconda, than all the dependencies are installed automatically except *cvxopt*. To install it see Dependencies
   If you dont have Anaconda installed and running ```python install.py``` failed, than the easy solution is to install Anaconda. Otherwise, search on the web how to install each dependency in the list appears in Dependencies.

### Citing GLINT
elior..todo


### Authors

This software was developed by Reut Yedidim, Omer Weissbrod  and Elior Rahmani.

For any question and for reporting bugs please send an email to Elior Rahmani at: elior.rahmani@gmail.com
