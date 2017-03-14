import os
import logging
from logging import handlers
import sys
import time
import traceback

LOG_FILE_PATH = 'glint.log' # TODO move to a configuration file
LOG_FILE_MAX_BYTES = 10 * 1024
LOG_FILE_BACKUP_COUNT = 5

OPTIONAL_LEVELS =   {"warning": logging.WARNING,
                     "critical":logging.CRITICAL,
                     "fatal":   logging.FATAL,
                     "info":    logging.INFO,
                     "debug":   logging.DEBUG,
                     "error":   logging.ERROR}


class _Formatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        super(_Formatter, self).__init__(*args, **kwargs)

    def get_record_info(self, record):
        source = '{}:{}'.format(record.filename, record.lineno)
        message = record.getMessage()
        timestamp = self.formatTime(record)
        exception = ''
        if record.exc_info is not None:
            exception_type, exception, traceback_obj = record.exc_info
            traceback_str = ' ---> '.join(['%s:%d, in %s: %s' % (filename, line_number, func_name, text) for (filename, line_number, func_name, text) in traceback.extract_tb(traceback_obj)])
            exception = '; EXCEPTION: %s, %s; TRACEBACK: %s' % (exception_type, exception, traceback_str)
        return source, message, timestamp, exception

class _FileFormatter(_Formatter):
    def format(self, record):
        source, message, timestamp, exception = super(_FileFormatter, self).get_record_info(record)
        return '{} glint {:<9} {:<30} {}{}'.format(timestamp, record.levelname, source, message, exception)


class _ConsoleFormatter(_Formatter):
    def format(self, record):
        source, message, timestamp, exception = super(_ConsoleFormatter, self).get_record_info(record)
        return '{:<9} {}{}'.format(record.levelname, message, exception)



class Configure(object):
    def __init__(self, loglevel=logging.INFO, prefix = ''):
        logging.raiseExceptions = 0
        logging.captureWarnings(True)

        self.logger = logging.getLogger()

        self.streamHandler = logging.StreamHandler()
        self.streamHandler.setFormatter(_ConsoleFormatter())
        self.logger.addHandler(self.streamHandler)

        self.setLoggerLevel(loglevel)
        if prefix:
            self.setLoggerFile(prefix)            

    def setLoggerLevel(self, loglevel):
        self.loglevel = loglevel
        self.logger.setLevel(loglevel)
        self.streamHandler.setLevel(loglevel)
    
    def setLoggerFile(self, prefix):
        filename = LOG_FILE_PATH
        if prefix:
            dirname = os.path.dirname(filename)
            basefilename = prefix + "." + os.path.basename(filename)
            filename = os.path.join(dirname, basefilename)
    
        self.fileHandler = logging.FileHandler(filename, mode='w')
        self.fileHandler.setFormatter(_FileFormatter())
        self.logger.addHandler(self.fileHandler)
        self.fileHandler.setLevel(self.loglevel)
