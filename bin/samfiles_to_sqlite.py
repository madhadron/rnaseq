"""
samfiles_to_sqlite.py
by Fred Ross, <madhadron@gmail.com>

Load a list of SAM/BAM files into an SQLite database as a group,
calculating their leftsite depth on each transcript and the
multiplicities.  The resulting database is ready to be handed to
find_subproblems.py and to inference.py, both of which will treat it
as read only.
"""

import getopt
import os
import sys
import pysam
import sqlite3
from rnaseq.load import *

usage = """samfiles_to_sqlite.py [-vh] [-l readlen] (-c|-x) [-g group] db samfiles ...

-v           Run verbosely
-h           Print this message and exit
-l readlen   Reads have length 'readlen' in SAM/BAM files
-c|-x        -c makes this group a control, -x makes it an experimental sample
-g group     Insert the samfiles into the database with group ID 'group'.  If
             omitted, just uses the next free value in the database.
db           The SQLite3 database to write to.
samfiles     A list of SAM/BAM files, each containing the aligned reads
             from one biological sample.
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class State(object):
    def __init__(self):
        self.verbose = False
        self.read_length = 38
        self.is_control = None
        self.group_id = None
        self.group_label = ""

state = State()

def vmsg(msg):
    if state.verbose:
        print >>sys.stderr, msg

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "hvl:L:cxg:", 
                                       ["help","read-length",
                                        "group-label","control",
                                        "experimental",
                                        "group-id"])
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
            elif o in ("-l", "--read-length"):
                try:
                    state.read_length = int(a)
                except ValueError, v:
                    raise Usage("Read length must be an integer, found %s" % a)
            elif o in ("-x", "--experimental"):
                state.is_control = False
            elif o in ("-c", "--control"):
                if state.is_control == False:
                    raise Usage("Can only specify one of -x and -c")
                else:
                    state.is_control = True
            elif o in ("-g", "--group-id"):
                try:
                    state.group_id = int(a)
                except ValueError, v:
                    raise Usage("Group ID must be an integer, found %s" % a)
            elif o in ("-L", "--group-label"):
                state.group_label = a
            else:
                raise Usage("Unhandled option: " + o)
        if len(args) < 2:
            raise Usage("samfiles_to_sqlite.py takes at least two arguments.")
        if state.is_control == None:
            raise Usage("Must specify one of -c or -x")

        db_filename = args[0]
        samfiles = args[1:]
        for f in samfiles:
            if not(os.path.exists(f)):
                raise Usage("Input file %s does not exist." % f)

        db_exists = os.path.exists(db_filename)
        db = sqlite3.connect(db_filename)
        if not(db_exists):
            initialize_database(db)

        control_sample_group = insert_sample_group(db, state.group_label, 
                                                   state.is_control, 
                                                   state.group_id)

        for f in samfiles:
            load_sam(db, f, control_sample_group)

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

