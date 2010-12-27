"""
find_subproblems.py
by Fred Ross, <madhadron@gmail.com>

Given an SQLite3 database as produced by samfiles_to_sqlite.py, find
groups of transcripts which can be treated as separate subproblems.
Each group of transcripts is printed as a space separated list of
integers referring to transcript IDs in the database.

This script uses a graph which is a union of all comparisons between
sample groups.  Biologically, this is unlikely to make much difference
since two transcripts which are connected in one pair of groups are
likely to be connected in most groups.
"""

import getopt
import os
import sys
import sqlite3
import rnaseq.subproblems

usage = """find_subproblems.py [-vh] db

-v   Run verbosely
-h   Print this message and exit
db   Database file to find subproblems in
"""

class Usage(Exception):
    def __init__(self,  msg):
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
        if len(args) != 1:
            raise Usage("find_subproblems.py takes exactly one argument.")

        db_filename = args[0]
        if not(os.path.exists(db_filename)):
            raise Usage("Database file %s does not exist" % db_filename)

        db = sqlite3.connect(db_filename)
        for q in find_subproblems(db):
            print ' '.join([str(x) for x in q])

        db.close()
    
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

