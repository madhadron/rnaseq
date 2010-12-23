"""
inference.py
by Fred Ross <madhadron@gmail.com>

Given a database of leftsites and multiplicities, and a list of
transcript IDs in that database forming an independent subproblem, do
one way linear model MCMC on two sample groups, pickle the resulting
posteriors, and print the names of pickle files.
"""

import getopt
import os
import sys
import pysam
import sqlite3
import pickle
from pymc import *
from numpy import *

import rnaseq

usage = """inferece.py [-vh] [-n nsamples] [-o picklefile] db group1 group2 1 2 3 ...

-v           Run verbosely
-h           Print this message and exit
nsamples     The number of MCMC samples to simulate
picklefile   The file to write posteriors to (defaults to db-group1-group2-1,2,3.pickle)
db           The database to read from
group1/2     Sample group IDs to use as the two conditions.
1 2 3 ...    A list of integers giving the transcripts to do inference on.
"""


class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class State(object):
    def __init__(self):
        self.verbose = False
        self.read_length = 38
        self.n_samples = 500
        self.pickle_filename = None

state = State()

def vmsg(msg):
    if state.verbose:
        print >>sys.stderr, msg

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "hvn:o:", ["help"])
        except getopt.error, message:
            raise Usage(message)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-v", ):
                state.verbose=True
                print "Running verbosely."
            elif o in ("-n", ):
                try:
                    state.n_samples = int(a)
                except ValueError, v:
                    raise Usage("Number of samples must be an integer, found %s" % a)
            elif o in ("-o",):
                state.pickle_filename = a
            else:
                raise Usage("Unhandled option: " + o)

        if len(args) < 4:
            raise Usage("inference.py takes at least four arguments.")
        db_filename = args[0]
        if not(os.path.exists(db_filename)):
            raise Usage("File %s does not exist" % db_filename)
        db = sqlite3.connect(db_filename)

        try:
            group1 = int(args[1])
        except ValueError, v:
            raise Usage("ID for group1 must be an integer, got %s" % args[1])
        
        try:
            group2 = int(args[2])
        except ValueError, v:
            raise Usage("ID for group2 must be an integer, got %s" % args[2])
        try:
            transcripts = [int(s) for s in args[3:]]
        except ValueError, v:
            raise Usage("Transcript IDs must be integers.")

        if state.pickle_filename == None:
            state.pickle_filename = \
                "%s-%d-%d-%s.pickle" % (db_filename, 
                                        group1, group2, 
                                        ','.join([str(i) for i in transcripts]))
        if os.path.exists(state.pickle_filename):
            raise Usage("Output file %s already exists.  Not overwriting." % state.pickle_filename)

        M = rnaseq.build_model(db, group1, group2, transcripts)

        M.sample(state.n_samples*5 + 2000, burn=2000, thin=5)

        with open(state.pickle_filename, 'w') as pf:
            [minusmu,a] = [{},{}]
            for i in transcripts:
                minusmu[i] = M.trace('minusmu'+str(i))[:]
                a[i] = M.trace('a'+str(i))[:]
            pickle.dump({'minusmu':minusmu, 'a':a}, pf)
    
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())
