# GLINT

GLINT is a user-friendly command-line tool for fast analysis of genome-wide DNA methylation data generated using the Illumina human methylation arrays. GLINT allows to easily run a pipeline of Epigenome-Wide Association Study (EWAS) under different models while accounting for known confounders in methylation data.

For more details about GLINT see the <a href="http://glint-epigenetics.readthedocs.io/" target="_blank">documentation</a>.

### Download and Installation

1. Download the latest release of GLINT from <a href="https://github.com/cozygene/glint/releases" target="_blank">here</a>.
2. Install the latest release of <a href="https://www.continuum.io/downloads" target="_blank">Anaconda for Python 2.7</a> which includes most of the dependencies required by GLINT. Note that GLINT requires Anaconda version >=4 (i.e. Anaconda2-4.X or above).
    - If you already have Python 2.7 and do not want to install Anaconda please see "Dependencies" bellow.
3. Use your terminal to install the *cvxopt* package for Python using Anaconda   
    **Windows** run ```conda install -c omnia cvxop```  
    **Linux or MacOS** run ```sudo `which conda` install -c anaconda cvxopt```   
    
### Documentation and a quick start tutorial
A detailed documentation of GLINT can be found <a href="http://glint-epigenetics.readthedocs.io/" target="_blank">here</a>. In addition, we provide a <a href="http://glint-epigenetics.readthedocs.io/en/latest/tutorial.html" target="_blank">quick start tutorial</a> that will get you started with GLINT quickly.

 
### Dependencies

GLINT was implemented for Python 2.7 and has the following dependencies:

    numpy (>=1.10.4)
    scipy (>=0.17)
    sklearn (>=0.17.1)
    pandas (>=0.18)
    matplotlib (>=1.5.1)
    statsmodels (>=0.6.1)
    cvxopt (>=1.1.8; not included in Anaconda by default)
    
We recommend installing the latest release of <a href="https://www.continuum.io/downloads" target="_blank">Anaconda for Python 2.7</a>, which already includes most of the necessary dependencies. In case you already have Python 2.7 installed and do not want to install Anaconda Python, run the "install.py" script we provide:
```
python install.py
```
The script automatically tries to install missing dependencies that are required for GLINT. Note that in some environments the script may fail to install some of the dependencies, in which case you will need to either manually install them or download <a href="https://www.continuum.io/downloads" target="_blank">Anaconda for Python 2.7</a>.


### Troubleshooting

#### Missing dependency?
1. Make sure you have Anaconda isntalled, see "Download and installation" for more details.  
  If you do not have Anaconda installed and running ```python install.py``` failed, then the easy solution is to install Anaconda. Otherwise, you will need to install all the Python packages listed under "Dependencies".
  If you have Anaconda installed, make sure you have updated versions of the required dependencies (as detailed above under "Denpendencies").
2. Make sure you run GLINT with the Anaconda Python command line (assuming you installed Anaconda Python):  
  **Windows**:  
    a. Find the path to Anaconda Python: press "Start" (win-key) and search for "conda". When you find it, do not open it but right click on it -> "properties", and extract the path written under "Location".  
    b. Use the path to the Anaconda Python instead of the standard "python" command. For example, if you found the path to be "C:\Users\me\Anaconda2\Scripts" then Anaconda Python is supposed to be under C:\Users\me\Anaconda2\python, in which case run each GLINT command as follows: ```C:\Users\me\Anaconda2\python glint.py...```  
  
  **Linux or MacOS**:  
    a. Find the path to Anaconda Python: run in the terminal: ```which conda```. If that command does not return anything then you do not have Anaconda installed; refer to "Download and Installation".  
    b. Use the path to the Anaconda Python instead of the standard "python" command. For example, if you found the path to be "/home/user/anaconda2/bin/" then Anaconda Python is supposed to be under /home/user/anaconda2/bin/python, in which case run each GLINT command as follows: ```/home/user/anaconda2/bin/python glint.py ...```  
 
#### MacOs cannot show .eps figures?
Restart the Preview app - quit it and try to open the figure again.

### Citing GLINT
If you use GLINT in any published work please cite it. Information about how to cite GLINT can be found in the documentation <a href="http://glint-epigenetics.readthedocs.io/en/latest/howtocite.html" target="_blank">here</a>.

### Authors

This software was developed by Reut Yedidim and Elior Rahmani. Code contributions were made by Omer Weissbrod and Dan coster. For any question and for reporting bugs or suggesting new features please send an email to Elior Rahmani at: elior.rahmani@gmail.com

