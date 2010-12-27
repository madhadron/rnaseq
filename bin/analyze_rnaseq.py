"""
analyze_rnaseq.py
by Fred Ross, <madhadron@gmail.com>

Load a configuration file specifying groups of FASTQ files, and
analyze the FASTQ files as RNASeq reads in pairs.
"""

import getopt
import sys
import os
from ConfigParser import SafeConfigParser, NoOptionError
from string import lstrip, rstrip
from bein import *

usage = """analyze.rnaseq.py [-hv] bowtie_index config_file

-h            Print this message and exit
-v            Run verbosely
bowtie_index  Index to run bowtie on the FASTQ files against
config_file   INI file specifying groups of FASTQ files
"""

def to_absolute(filename):
    abs_filename = os.path.abspath(filename)
    if os.path.exists(abs_filename):
        return abs_filename
    else:
        raise ConfigurationError("File %s in configuration doesn't exist." % filename)

def load_configuration(config_file):
    cp = SafeConfigParser({})
    cp.read(config_file)
    vmsg("Loaded configuration from %s" % config_file)
    r = {}
    for s in cp.sections():
        try:
            files = [lstrip(rstrip(f)) for f in 
                     cp.get(s, 'fastqfiles').split(',')]
            abs_files = [to_absolute(f) for f in files]
            r[s] = {'control': cp.getboolean(s, 'control'),
                    'fastqfiles': abs_files}
        except NoOptionError, n:
            raise ConfigurationError("Group '%s' in %s is missing '%s' option" % (n.section,
                                                                                  config_file,
                                                                                  n.option))
        vmsg("Added group %s" % s)
    return r

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class ConfigurationError(Exception):
    def __init__(self, msg):
        self.msg = msg

class State(object):
    def __init__(self):
        self.verbose = False

state = State()

def vmsg(msg):
    if state.verbose:
        print >>sys.stderr, msg

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "hv", 
                                       ["help","verbose"])
        except getopt.error, message:
            raise Usage(message)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-v", "--verbose"):
                state.verbose=True
                print "Running verbosely."
            else:
                raise Usage("Unhandled option: " + o)
        if len(args) != 2:
            raise Usage("analyze_rnaseq.py takes exactly two arguments.")

        bowtie_index = args[0]
        config_file = args[1]
        if not(os.path.exists(config_file)):
            raise Usage("Configuration file %s does not exist." % config_file)

        groups = load_configuration(config_file)


        with execution(None) as ex:
            # Make parallel bowtie recurse on dictionaries and nested lists, and return values in the same structure
        # in execution:
        # * run bowtie in parallel on fastq files
        # * Load sam files into a database
        # * Find subproblems
        # * Calculate pairs of groups to use
        # * Do pairwise inference on all subproblems
        # * Merge all pickles into results database and add to repository
        # * Calculate means and 95% CIs for all groups/genes and write as text file

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    except ConfigurationError, cfe:
        print >>sys.stderr, "Error in configuration file:"
        print >>sys.stderr, cfe.msg
        return 2

if __name__ == '__main__':
    sys.exit(main())

