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
    def __init__(self, namespace, *args, **kwargs):
        self._namespace = namespace
        super(_Formatter, self).__init__(*args, **kwargs)

    def formatException(self, record):
        if record.exc_info is not None:
            exceptionType, exception, tracebackObject = record.exc_info
            tracebackString = ' ---> '.join(['%s:%d, in %s: %s' % (filename, lineNumber, functionName, text) for (filename, lineNumber, functionName, text) in traceback.extract_tb(tracebackObject)])
            return '; EXCEPTION: %s, %s; TRACEBACK: %s' % (exceptionType, exception, tracebackString)
        else:
            return ''

class _FileFormatter(_Formatter):
    def format(self, record):
        timestamp = self.formatTime(record)
        source = '{}:{}'.format(record.filename, record.lineno)
        message = record.getMessage()
        exception = self.formatException(record)
        
        return '{} glint {:<3} {:<8} {:<32} {}{}'.format(timestamp, self._namespace, record.levelname, source, message, exception)


class _ConsoleFormatter(_Formatter):
    def format(self, record):
        source = '{}:{}'.format(record.filename, record.lineno)
        message = record.getMessage()
        exception = self.formatException(record)
        
        return '{:<10} {}{}'.format(record.levelname, message, exception)



class ConfigureLogging(object):
    def __init__(self, loglevel=logging.INFO, namespace='', prefix = ''):
        logging.raiseExceptions = 0
        logging.captureWarnings(True)

        self.logger = logging.getLogger()

        self.streamHandler = logging.StreamHandler()
        self.streamHandler.setFormatter(_ConsoleFormatter(namespace))
        self.logger.addHandler(self.streamHandler)


        self.fileHandler = logging.FileHandler(os.path.join(os.path.dirname(LOG_FILE_PATH), prefix + os.path.basename(LOG_FILE_PATH)),
                                          mode='w')
        self.fileHandler.setFormatter(_FileFormatter(namespace))
        self.logger.addHandler(self.fileHandler)

        self.setLoggerLevel(loglevel)

    def setLoggerLevel(self, loglevel):
        self.logger.setLevel(loglevel)
        self.streamHandler.setLevel(loglevel)
        self.fileHandler.setLevel(loglevel)
