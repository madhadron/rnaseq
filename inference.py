"""
inference.py
by Fred Ross <madhadron@gmail.com>

Given a database of leftsites and multiplicities, and a list of
transcript IDs in that database forming an independent subproblem,
doing one way linear model MCMC, pickle the resulting posteriors, and
print the names of pickle files.
"""

import getopt
import os
import sys
import pysam
import sqlite3
from pymc import *
from numpy import *

import rnaseq

usage = """inferece.py [-vh] db 1 2 3 ...

-v           Run verbosely
-h           Print this message and exit
db           The database to read from
1 2 3 ...    A list of integers giving the transcripts to do inference on.
"""


class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class State(object):
    def __init__(self):
        self.verbose = False
        self.read_length = 38

state = State()

def vmsg(msg):
    if state.verbose:
        print >>sys.stderr, msg

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "hv", ["help"])
        except getopt.error, message:
            raise Usage(message)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            if o in ("-v", ):
                state.verbose=True
                print "Running verbosely."
            else:
                raise Usage("Unhandled option: " + o)

        if len(args) < 2:
            raise Usage("inference.py takes at least two arguments.")
        db_filename = args[0]
        if not(os.path.exists(db_filename)):
            raise Usage("File %s does not exist" % db_filename)
        db = sqlite3.connect(db_filename)

        try:
            transcripts = [int(s) for s in args[1:]]
        except ValueError, v:
            raise Usage("Transcript IDs must be integers.")
        M = rnaseq.build_model(db, transcripts)
    
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
