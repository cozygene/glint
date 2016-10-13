# GLINT

GLINT is a user-friendly command line tool for fast analysis of genome-wide DNA methylation data generated using the Illumina human methylation arrays. GLINT allows to easily run a pipeline of Epigenome-Wide Association Study (EWAS) under different models while accounting for known confounders in methylation data.

For more details about GLINT see the <a href="blank" target="_blank">documentation</a>.

### Download and Installation

1. Download the latest release of GLINT from <a href="https://github.com/cozygene/glint/releases" target="_blank">here</a>.
2. Install the latest release of <a href="https://www.continuum.io/downloads" target="_blank">Anaconda for Python 2.7</a> which includes most of our dependencies.  
    - If you already have Python 2.7 and do not want to install Anaconda please see "Dependencies" bellow.
3. Install the *cvxopt* package for Python using Anaconda:   
    **Windows** run ```conda install -c omnia cvxop```  
    **Linux and MacOS** run ```sudo `which conda` install -c anaconda cvxopt```   
    
### Documentation and a quick start tutorial
A detailed documentation of glint can be found <a href="todo add link to docs" target="_blank">here</a>. In addition, we provide a <a href="todo add link to tutorial" target="_blank">quick start tutorial</a> that will get you started with GLINT quickly.

 
### Dependencies

GLINT was implemented for Python 2.7 and has the following dependencies:

    numpy
    scipy
    sklearn
    pandas
    matplotlib
    statsmodels
    cvxopt (not included in Anaconda by default)
    

We recommend installing the latest release of <a href="https://www.continuum.io/downloads" target="_blank">Anaconda for Python 2.7</a>, which already includes most of necessary dependencies. 
In order to install the dependency *cvxopt* with Anaconda run:

**Windows** : ```conda install -c omnia cvxop```

**Linux**: ```sudo `which conda` install -c anaconda cvxopt```

**MacOS**: ```todo```

In case the installation of *cvxopt* fails, we advise to search online: "how to install cvxopt with anaconda on [one of Windows/Linux/MacOS]"
 
In case you already have Python 2.7 installed and do not want to install Anaconda Python, run the "install.py" script we provide:
```
python install.py
```
The script automatically tries to install missing dependencies that are required for GLINT. Note that in some environments the script may fail to install some of the dependencies, in which case you will need to either manually install them or download <a href="https://www.continuum.io/downloads" target="_blank">Anaconda for Python 2.7</a>.

### Troubleshooting

#### Missing Dependency?
1. Make sure you have Anaconda isntalled, see Download and installation for more details.  
  If you dont have Anaconda installed and running ```python install.py``` failed, than the easy solution is to install Anaconda. Otherwise, search on the web how to install each dependency in the list appears in Dependencies.
2. Make sure you run GLINT with Anaconda Python command line:  
  **Windows**:  
    a. Find the path to python command line tool: press "Start" (win-key) and search for "conda". When you find it, dont open it but right click on it -> "properties" and there you can see the path under "Location".
   For this example, lets assume you found it at "C:\Users\me\Anaconda2\Scripts" so Anaconda Python supposed to be at C:\Users\me\Anaconda2\python  
    b. Run GLINT with that path: run ```C:\Users\me\Anaconda2\python glint.py...```  
  
  **Linux or MacOS**: 
    a. Find the path to python command line tool: run on command line: ```which conda```. If that command returns nothing than you don't have Anaconda installed, refer to Download and Installation. Otherwise, if for example the output of the command is /home/user/anaconda2/bin/conda than Anaconda Python command line is at /home/user/anaconda2/bin/python.  
    b. run GLINT with that path: run ```/home/user/anaconda2/bin/python glint.py ...```  
 
#### MacOs can't show .eps Figure?
Restart the Preveiw app: quit it and try to open the figure again.

### Citing GLINT
If you use GLINT in any published work please cite it. Information about how to cite GLINT can be found in the documentation <a href="howtocite.html" target="_blank">here</a>.

### Authors

This software was developed by Reut Yedidim and Elior Rahmani. Code contributions were made by Omer Weissbrod and Dan coster. For any question and for reporting bugs or suggesting new features please send an email to Elior Rahmani at: elior.rahmani@gmail.com
