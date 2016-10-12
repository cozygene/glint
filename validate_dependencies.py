import sys
import os
from distutils import spawn
from install import check_dependencies
import subprocess

"""
doc:
this script tries to import all dependencies packages. if it fails there are two options:
1 - user didn't run install.py (the script which installs dependencies)
2 - user installed anaconda after he already had python installed (shouldn't happen on linux )

this script:
1 - tries to import dependencies. on succes - exit.
2 - if it failes, on linux - tells the user to run installation script. on windows it searches for anaconda - if it doesnt find it - it tells the user to run install.py (which will tell him to install anaconda)
3 - if it fount anaconda, it tries to add it to the path and then import dependencies - if it fails it tells the user how to add it to the path by himself
"""

# dependencies that included in anaconda
GLINT_OBLIGATORY_DEPENDENCIES_WITH_CONDA = ['numpy', 'scipy', 'sklearn', 'matplotlib', 'pandas', 'statsmodels'] # TODO move to configuration file
# dependencies that are not included in anaconda
GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA = ['cvxopt'] # TODO move to configuration file
INSTALL_WINDODWS = {'cvxopt':"echo hi"}
INSTALL_LINUX = {'cvxopt': "sudo `which conda` install -c anaconda cvxopt"} #conda install -c anaconda cvxopt=1.1.8

ERROR_MSG = "Something is wrong with GLINT dependencies, please follow the README at https://github.com/cozygene/glint"

PYTHONPATH_EXPLAIN = \
"""
Something is wrong with the environment. This probably happened since you had python installed on your pc before you installed anaconda.
All you need to do is add the Anaconda path to your environment variables. Follow these steps:

1. find the path where Anaconda was installed. One way to do that is to press "Start" (win-key) and search for "conda",
   when you find it, dont open it but right click on it -> "properties" and there you can see the path.
   For this example, lets assume you found it at "C:\Users\me\Anaconda".

2. look for "site-packages"  folder insite Anaconda folder you found (you can use the "search" box in the explorer)
    For this example, lets assume you found it at "C:\Users\me\Anaconda\Lib\site-packages"

3. add the full "site-packages" path to system variables:
    3.1 go to My Computer > Properties > Advanced System Settings > Environment Variables 
    3.2 edit PYTHONPATH variable: insert the "site-packages" path at the begining
        for this example, lets assume PYTHONPATH is now "C:\Python27\Lib;C:\Python27\DLLs;"
        change it to "C:\Users\me\Anaconda\Lib\site-packages;C:\Python27\Lib;C:\Python27\DLLs;"
4. try to run glint again
"""
DEPENDENCIES_EXPLANATION = "\nsome dependencies are missing.\nto install either run\n\tpython install.py\nor install anaconda from https://www.continuum.io/downloads.\ndependencies are: %s\nplease contact us on failure"
INSTALL_ANACONDA = "\nSome dependencies are missing.\nTo install them please install Anaconda Pythn from https://www.continuum.io/downloads.\ndependencies are: %s\nplease contact us on failure\nFor more information see the README at https://github.com/cozygene/glint" 
RUN_INSTALL_SCRIPT = "please run install.py"

def import_dependencies(dependencies_list):
    for module_name in dependencies_list:
        __import__(module_name)
    
def install_without_anaconda():
    if sys.platform.startswith("win") or "nt" in os.name: # you run windows, we dont know to install dependencie sthere
        print INSTALL_ANACONDA
        sys.exit()
    elif sys.platform.startswith("lin") or "posix" in os.name: # you run linux
        # install dependencis that are not installed
        not_installed = check_dependencies(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA + GLINT_OBLIGATORY_DEPENDENCIES_WITH_CONDA)
        if not_installed == []: # we succeeded to install all dependencies
            return True
        print INSTALL_ANACONDA
        sys.exit()

def install_packages_with_conda(dependencies_list):
    install_instructions = dict()
    if sys.platform.startswith("win") or "nt" in os.name: # you run windows
        install_instructions = INSTALL_WINDODWS
    elif sys.platform.startswith("lin") or "posix" in os.name: # you run linux
        install_instructions = INSTALL_LINUX
    
    for pkg in dependencies_list:
        if pkg in install_instructions:
            res = subprocess.check_output(install_instructions[pkg], shell=True, stderr = subprocess.STDOUT)
            try:
                import_dependencies([pkg])
            except:
                print "Could not install package %s" % pkg
                return False
        else:
            return False
    return True

def add_anaconda_to_path_lin(path):
    """
    func assumes anaconda installed on linux
    """
    anaconda_dir = os.path.dirname(path).lower()
    conda_python = os.path.join(anaconda_dir, "python")
    if os.path.exists(conda_python):
        print "You have Anaconda installed, switching to Anaconda..."
        os.execv(conda_python, [conda_python] + sys.argv)

def add_anaconda_to_path_win(path):
    """
    only for windows
    search for anaconda on the PC
    try to add it's path to PYTHONPATH and import dependencies again
    if it failes - it instruct the user what to do next
    """

    # if anaconda was found, try to add it to path and import dependencies
    
    success = False
    if os.path.exists(os.path.join(path,"Lib","site-packages")):
        module_dir = os.path.join(path,"Lib","site-packages")
        try:
            sys.path.insert(0, module_dir) # add to PYTHONPATH 
            import_dependencies()
            success = True
        except:
            del sys.path[0]

    # if we couldn't add Anaconda to PATH, tell the user how to do it
    if not success:
        print PYTHONPATH_EXPLAIN

def run_me_with_anaconda():
    if sys.platform.startswith("win") or "nt" in os.name: # you run windows
        path = os.path.dirname(conda_path)
        last_path = ""
        while last_path != path and "conda" not in os.path.basename(path).lower():
            last_path = path
            path = os.path.dirname(path)
        add_anaconda_to_path_win(path)
    elif sys.platform.startswith("lin") or "posix" in os.name: # you run linux
        add_anaconda_to_path_lin(conda_path)
    else:
        print "Please read the readme at https://github.com/cozygene/glint"
        sys.exit()
        # DEPENDENCIES_EXPLANATION % "\n\t".join(GLINT_OBLIGATORY_DEPENDENCIES)

print "Validating all dependencies are installed..."
# search if anaconda exist
conda_path = spawn.find_executable("conda") # "conda" is the command line installed by anaconda
# if anaconda wasnt found, tell the user to run installation script
if not conda_path:
    install_without_anaconda()

# see if we are running anaconda python
anaconda_dir = os.path.dirname(conda_path).lower()
if anaconda_dir in sys.executable.lower():
    print "You are now running Anaconda Python"
    try:
        for pkg in GLINT_OBLIGATORY_DEPENDENCIES_WITH_CONDA:
            import_dependencies([pkg])
    except:
        # something is wrong since we are missing dependencies that included in anaconda
        print ERROR_MSG
        print pkg
        sys.exit()

    try:
        import_dependencies(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA)
    except:
        print "Some dependencies are missing, trying to install them..."
        success = install_packages_with_conda(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA)
        if not success:
            print ERROR_MSG
            sys.exit()
else: 
    # we are not running anaconda python - execute it
    run_me_with_anaconda()

print "All dependencies are installed"
