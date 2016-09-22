import sys
import os
from distutils import spawn

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

GLINT_OBLIGATORY_DEPENDENCIES = ['numpy', 'scipy', 'sklearn', 'matplotlib', 'pandas', 'statsmodels'] # TODO move to configuration file

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

RUN_INSTALL_SCRIPT = "please run install.py"
def import_dependencies():
    for module_name in GLINT_OBLIGATORY_DEPENDENCIES:
        __import__(module_name)
    
def add_anaconda_to_path():
    """
    only for windows
    search for anaconda on the PC
    try to add it's path to PYTHONPATH and import dependencies again
    if it failes - it instruct the user what to do next
    """
    conda_path = spawn.find_executable("conda") # "conda" is the command line installed by anaconda
    
    # if anaconda wasnt found, tell the user to run installation script
    if not conda_path:
        print RUN_INSTALL_SCRIPT
        return 

    # if anaconda was found, try to add it to path and import dependencies
    path = os.path.dirname(conda_path)
    last_path = ""
    while last_path != path and "conda" not in os.path.basename(path).lower():
        last_path = path
        path = os.path.dirname(path)
    
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


    # try to import dependencies
    # if it fails, search for anaconda and 
try:
    import_dependencies()
except:
    if "win" in sys.platform or "nt" in os.name: # you run windows
        add_anaconda_to_path()
    else:
        print 3
        print RUN_INSTALL_SCRIPT
        print "if you just installed anaconda, try to open new command line and run glint again. contact us on failure."
    sys.exit()
