#!python
"""
inference_subproblem.py
by Fred Ross, <madhadron@gmail.com>

Runs a one-way linear model on a given set of transcripts as loaded in a given database.  Writes the result as a pikle of a dictionary of parameters pointing to a dictionary of transcript keys pointing to NumPy arrays of the posteriors for those samples.
"""

import getopt
import os
import sys
import sqlite3
import pickle
from rnaseq import *

usage = """inference_subproblem.py [-vh] [-n n_samples] pickle_file db group1 group2 transcripts ...

-v             Run verbosely
-h             Print this message and exit
db             The SQLite3 database to read from.
-n n_samples   Produce n_samples samples of the posterior.
group1,group2  Integers giving the group IDs to work on.
transcripts    Integers giving the transcripts to do inference on.
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class State(object):
    def __init__(self):
        self.verbose = False
        self.n_samples = 500

state = State()

def vmsg(msg):
    if state.verbose:
        print >>sys.stderr, msg

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "hvn:", ["help","verbose"])
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
            elif o in ("-n", ):
                try:
                    state.n_samples = int(a)
                except ValueError, v:
                    raise Usage("Number of samples must be an integer, found %s" % a)
            else:
                raise Usage("Unhandled option: " + o)
        if len(args) < 5:
            raise Usage("simple_inference.py takes at least five arguments.")

        pickle_filename = args[0]
        if os.path.exists(pickle_filename):
            raise Usage("Pickle file %s already exists." % pickle_filename)

        db_filename = args[1]
        if not(os.path.exists(db_filename)):
            raise Usage("No such database %s" % db_filename)
        else:
            db = sqlite3.connect(db_filename)

        try:
            group1 = int(args[2])
        except ValueError, v:
            raise Usage("group1 must be an integer; found %s" % args[2])

        try:
            group2 = int(args[3])
        except ValueError, v:
            raise Usage("group2 must be an integer; found %s" % args[3])

        try:
            transcripts = [int(x) for x in args[4:]]
        except ValueError, v:
            raise Usage("All transcripts must be integers; found %s" %
                        ', '.join([str(t) for t in args[4:]]))

        # Check that the given transcripts form a complete subproblem
        if len(transcripts) > 1:
            query = """select distinct a.transcript
                       from (select * from multiplicity_entries
                             where transcript not in %s) as a
                       join (select * from multiplicity_entries 
                             where transcript in %s) as b
                       on a.multiplicity = b.multiplicity
                    """ % (str(tuple(transcripts)), str(tuple(transcripts)))
            missed_transcripts = [x for (x,) in db.execute(query)]
            if missed_transcripts != []:
                raise Usage("The given set of transcripts, %s, is incomplete.  The complete subproblem also contains %s" % (', '.join([str(t) for t in transcripts]),
                                                                                                                            ', '.join([str(t) for t in missed_transcripts])))
                      

        vmsg("Doing inference on transcripts %s" %
             ', '.join([str(t) for t in transcripts]))
        M = build_model(db, group1, group2, transcripts)
        vmsg("Built model")
        M.sample(state.n_samples*5 + 2000, burn=2000, thin=5)
        vmsg("Sampled from model")

        posteriors = {}
        for i in transcripts:
            posteriors[i] = {'mu': -1*M.trace('minusmu'+str(i))[:],
                             'a': M.trace('a'+str(i))[:]}
        with open(pickle_filename, 'w') as pf:
            pickle.dump((group1,group2,posteriors), pf)

        vmsg("Wrote pickle file in %s" % pickle_filename)

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

