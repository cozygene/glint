import logging
import os

CUR_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
TUTORIAL_FILES_DIR = "tutorial"

TUTORIAL_FILES = ["tutorial/datafile.txt", "tutorial/covariates.txt", "tutorial/phenotypes.txt"]

def check_tutorial_files_exist():
    for path in TUTORIAL_FILES:
        if not os.path.exists(path):
            print "make sure you have the file %s under the '%s' folder" % (os.path.basename(path),TUTORIAL_FILES_DIR)
            exit(2)

class TutorialTester():
    TUTORIAL_CMDS = [
                    "python glint.py --datafile tutorial/datafile.txt --covarfile tutorial/covariates.txt --phenofile tutorial/phenotypes.txt --gsave",
                    "python glint.py --datafile datafile.glint --plot --plotpcs --numpcs 2 --out pcs_plot",
                    "python glint.py --datafile datafile.glint --maxpcstd 1 4 --gsave --out data_cleaned",
                    "python glint.py --datafile data_cleaned.glint --refactor --k 6 --covar age gender chip1 chip2 chip3 chip4 chip5 chip6 chip7 chip8 --gsave --out data_cleaned_v2",
                    "python glint.py --datafile data_cleaned_v2.glint --epi --covar rc1 rc2 rc3 rc4 rc5 rc6 --gsave --out data_final",
                    "python glint.py --datafile data_final.glint --ewas --linreg --pheno y1 --covar age gender rc1 rc2 rc3 rc4 rc5 rc6 epi1 --stdth 0.01 --rmxy --rmns --rmpoly",
                    "python glint.py --plot --qqplot --manhattan --results results.glint.linreg.txt",
                    "python glint.py --datafile data_final.glint --ewas --linreg --pheno y1 --covar age gender epi1 --stdth 0.01 --rmxy --rmns --rmpoly --plot --qqplot --manhattan --out unadjusted",
                    "python glint.py --datafile tutorial/datafile.txt --gsave --out newdata",
                    "python glint.py --datafile datafile.glint --txtsave",
                    "python replace_missing_values.py --datafile tutorial/datafile.txt --chr NA --maxs 0.03 --maxi 0.03",
                    ]
    def __init__(self):
        logging.info("Testing Started on TutorialTester")
        check_tutorial_files_exist()
        self.run_tutorial_commands()
        logging.info("Testing Finished on TutorialTester")

    def run_tutorial_commands(self):
        for cmd in self.TUTORIAL_CMDS:
            logging.info("running: %s" , cmd)
            res = os.system(cmd)
            if res != 0:
                logging.error("error in tutorial command: '%s'" % cmd)
                exit(2)
        print "PASS"
