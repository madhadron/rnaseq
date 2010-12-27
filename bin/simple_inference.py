#!python
"""
simple_inference.py
by Fred Ross, <madhadron@gmail.com>

Script for running a one-way linear model on two sets of SAM/BAM files
of aligned reads.  The transcripts in each file must be identical.
The result is an SQLite database containing all the intermediate work,
plus the final posteriors from the analysis.
"""

import getopt
import os
import sys
import sqlite3
from rnaseq import *

usage = """simple_inference.py [-vh] db group1 group2

-v             Run verbosely
-h             Print this message and exit
db             The SQLite3 database to write to.
group1,group2  Comma separated list of SAM/BAM files to use as samples
               for the two conditions
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
            elif o in ("-v", ):
                state.verbose=True
                print "Running verbosely."
            elif o in ("-n", ):
                try:
                    state.n_samples = int(a)
                except ValueError, v:
                    raise Usage("Number of samples must be an integer, found %s" % a)
            else:
                raise Usage("Unhandled option: " + o)
        if len(args) < 3:
            raise Usage("simple_inference.py takes at least three arguments.")

        db_filename = args[0]
        if os.path.exists(db_filename):
            raise Usage("Database file %s already exists." % db_filename)
        db = sqlite3.connect(db_filename)
        initialize_database(db)
        vmsg("Initialized database %s" % db_filename)

        group1_files = args[1].split(',')
        group2_files = args[2].split(',')
        for f in group1_files + group2_files:
            if not(os.path.exists(f)):
                raise Usage("Input file %s does not exist." % f)

        group1_id = insert_sample_group(db, 'Group 1', False)
        group2_id = insert_sample_group(db, 'Group 2', False)

        vmsg("Adding group 1 files:")
        for f in group1_files:
            load_sam(db, f, group1_id)
            vmsg("Loaded %s" % f)

        vmsg("Adding group 2 files:")
        for f in group2_files:
            load_sam(db, f, group2_id)
            vmsg("Loaded %s" % f)

        db.execute("""insert into inferences (id,group1,group2) 
                      values (1,?,?)""",
                   (group1_id, group2_id))
        db.commit()

        subproblems = find_subproblems(db)

        for transcripts in subproblems:
            vmsg("Doing inference on transcripts %s" % ', '.join([str(t) for t in transcripts]))
            M = build_model(db, group1_id, group2_id, transcripts)
            M.sample(state.n_samples*5 + 2000, burn=2000, thin=5)
            for t in transcripts:
                mm = M.trace('minusmu%d' % t)[:]
                for i,v in enumerate(mm):
                    db.execute("""insert into posterior_samples
                                  (inference,transcript,variable,sample,value)
                                  values (1,?,'mu',?,?)""",
                               (t,i,-1*v))
                a = M.trace('a%d' % t)[:]
                for i,v in enumerate(a):
                    db.execute("""insert into posterior_samples
                                  (inference,transcript,variable,sample,value)
                                  values (1,?,'a',?,?)""",
                               (t,i,v))
                db.commit()

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

