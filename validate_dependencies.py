import sys
import os
from distutils import spawn
from install import check_dependencies, color_print, FOREGROUND, BACKGROUND
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
INSTALL_WINDODWS = {'cvxopt':"conda install -c omnia cvxop"}
INSTALL_LINUX = {'cvxopt': "sudo `which conda` install -c anaconda cvxopt"} #conda install -c anaconda cvxopt=1.1.8

DEPENDENCIES_EXPLANATION = "\nsome dependencies are missing.\nto install either run\n\tpython install.py\nor install anaconda from https://www.continuum.io/downloads.\ndependencies are: %s\nplease contact us on failure"
INSTALL_ANACONDA = "Some dependencies are missing"
TROUBLESHOOT = "Please follow the instructions at https://github.com/cozygene/glint/blob/master/README.md#troubleshooting"

def import_dependencies(dependencies_list):
    for module_name in dependencies_list:
        __import__(module_name)
    
def install_without_anaconda():
    if sys.platform.startswith("win") or "nt" in os.name: # you run windows, we dont know to install dependencie sthere
        color_print(INSTALL_ANACONDA +" "+ TROUBLESHOOT, FOREGROUND.RED)
        sys.exit()
    elif sys.platform.startswith("lin") or "posix" in os.name: # you run linux
        # install dependencis that are not installed
        not_installed = check_dependencies(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA + GLINT_OBLIGATORY_DEPENDENCIES_WITH_CONDA)
        if not_installed == []: # we succeeded to install all dependencies
            print("Dependencies are installed")
            return True
        color_print(INSTALL_ANACONDA+" "+ TROUBLESHOOT, FOREGROUND.RED)
        sys.exit()

def install_packages_with_conda(dependencies_list):
    install_instructions = dict()
    if sys.platform.startswith("win") or "nt" in os.name: # you run windows
        install_instructions = INSTALL_WINDODWS
    elif "posix" in os.name: # you run linux/macos 
        install_instructions = INSTALL_LINUX
    color_print("To install dependency '%s' please run:\n  %s" % (dependencies_list[0], install_instructions[dependencies_list[0]]), FOREGROUND.RED)
    sys.exit()

    # that code supposed
    # print("Trying to install dependencies...")
    # for pkg in dependencies_list:
    #     if pkg in install_instructions:
    #         try:
    #             res = subprocess.check_output(install_instructions[pkg], shell=True, stderr = subprocess.STDOUT)
    #         except Exception as e:
    #             print e
    #             color_print("Could not install cvxopt, " + TROUBLESHOOT, FOREGROUND.RED)
    #             sys.exit()
    #         try:
    #             import_dependencies([pkg])
    #         except:
    #             print "Could not install package %s" % pkg
    #             color_print("To install package '%s' please run:\n  %s" % (pkg, install_instructions[pkg]), FOREGROUND.RED)
    #             return False
    #     else:
    #         return False
    # return True

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
    anaconda_python = os.path.join(path,"python")
    color_print("\nYou have Anaconda installed but you are not running it's python, please run:\n\n%s %s" %(anaconda_python, " ".join(sys.argv)), FOREGROUND.RED)
    sys.exit()
    # success = False
    # if os.path.exists(os.path.join(path,"Lib","site-packages")):
        # module_dir = os.path.join(path,"Lib","site-packages")
        # sys.path.insert(0, r"C:\Users\me\Anaconda2\Lib\site-packages\statsmodels\base") # add to PYTHONPATH 

        # try:
            # import_dependencies(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA)
        # except:
			# print "Some dependencies are missing, trying to install them..."
			# success = install_packages_with_conda(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA)
			# if not success:
				# print ERROR_MSG, TROUBLESHOOT
				# sys.exit()
			
        # try:
            # import_dependencies(GLINT_OBLIGATORY_DEPENDENCIES_WITH_CONDA)
            # success = True
        # except:
            # del sys.path[0]
            # print "Some dependencies are missing,", TROUBLESHOOT
            # sys.exit()

def run_me_with_anaconda():
    if sys.platform.startswith("win") or "nt" in os.name: # you run windows
        path = os.path.dirname(conda_path)
        last_path = ""
        while last_path != path and "conda" not in os.path.basename(path).lower():
            last_path = path
            path = os.path.dirname(path)
        add_anaconda_to_path_win(path)
    elif "posix" in os.name: # you run linux/macos
        add_anaconda_to_path_lin(conda_path)
    else:
        color_print("You have Anaconda installed but you are not running it's python "+ TROUBLESHOOT, FOREGROUND.RED)
        sys.exit()

print "Validating all dependencies are installed..."
# search if anaconda exist
conda_path = spawn.find_executable("conda") # "conda" is the command line installed by anaconda
# if anaconda wasnt found, tell the user to run installation script
if not conda_path:
    install_without_anaconda()
    sys.exit()

#from now code assumes anaconda is installes
		
# see if we are running anaconda python
anaconda_dir = os.path.dirname(conda_path).lower()

if anaconda_dir in sys.executable.lower() or os.path.dirname(anaconda_dir) in sys.executable.lower():
    print "You are now running Anaconda Python"
    try:
        for pkg in GLINT_OBLIGATORY_DEPENDENCIES_WITH_CONDA:
            import_dependencies([pkg])
    except:
        # something is wrong since we are missing dependencies that included in anaconda
        color_print("There was a problem with the dependency %s, " % pkg + TROUBLESHOOT, FOREGROUND.RED)
        sys.exit()

    try:
        import_dependencies(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA)
    except:
        print "Some dependencies are missing"
        success = install_packages_with_conda(GLINT_OBLIGATORY_DEPENDENCIES_NO_CONDA)
        if not success:
            color_print("Could not install dependencies, ", TROUBLESHOOT, FOREGROUND.RED)
            sys.exit()
			
else: 
    # we are not running anaconda python - execute it
    run_me_with_anaconda()

print "All dependencies are installed"
