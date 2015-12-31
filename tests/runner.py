import logging
logging.basicConfig( format = '%(levelname)s:%(message)s', level = logging.DEBUG )
import subprocess
import os
import sys
import argparse 
import time
import tableprinter
TESTS_DIR = "."

class Runner( object ):
    def __init__( self, tests ):
        try:
            self._prepare()
            self._runTests( tests )
        except Exception:
            logging.exception( 'in tests runner' )
            raise

        self._reportTests()


    def _prepare( self ):
        self._success = True
        if not os.path.exists('logs'):
            os.mkdir('logs')
        self._onlineReport = open( 'logs/systemtest_%s.summary' % time.strftime( '%H_%M' ), 'w' )

    def _runTests( self, tests ):
        self._results = []
        for test in tests:
            logging.info( "*" * 80 )
            logging.info( "STARTING UP SYSTEM FOR %s" % test )
            logging.info( "RUNNING %s" % test )
            exitCode = subprocess.call( 'python %s' % test, shell = True )
            if exitCode != 0:
                exitString = 'FAILURE'
                self._success = False
            else:
                exitString = 'SUCCESS'
            logging.info( "RESULT OF %s -----------------------> %s" % ( test, exitString ) )
            self._onlineReport.write( '%s: %s\n' % ( exitString, test ) )
            self._onlineReport.flush()
            self._results.append( ( test, exitString ) )


    def _reportTests( self ):
        print '---- FINAL REPORT ----'
        head = [ 'test', 'result' ]
        rows = self._results
        tableprinter.TablePrinter( rows, head, 'll' ).output()
        if not self._success:
            print 'SYSTEM TEST FAILURE!'
            quit( 1 )
        else:
            print 'SYSTEM TEST SUCCESS!'

def find_all_test_files():
    tests = []
    for root, dirs, files in os.walk(TESTS_DIR):
        for file in files:
            if file.startswith("test_"):
                tests.append(os.path.join(root,file))
    return tests

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument( '--tests', nargs = '+' )
    arguments = parser.parse_args()

    if arguments.tests is None:
        tests = find_all_test_files()
    else:
        print "Running specified tests is NOT SUPPORTED right now"
    logging.info( "found tests: %s" % " ".join( tests ) )
    Runner( tests )
