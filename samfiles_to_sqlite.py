"""
samfiles_to_sqlite.py
by Fred Ross, <madhadron@gmail.com>

Load two lists of SAM/BAM files into an SQLite database, calculating
their leftsite depth on each transcript and the multiplicities.  The
resulting database is ready to be handed to find_subproblems.py and to
inference.py, both of which will treat it as read only.
"""

import getopt
import os
import sys
import pysam
import sqlite3

usage = """samfiles_to_sqlite.py [-vh] [-l readlen] condition1 condition2 db

-v           Run verbosely
-h           Print this message and exit
-l readlen   Reads have length 'readlen' in SAM/BAM files
condition1/2 Comma separated lists of SAM/BAM files, each containing
             the aligned reads from one biological sample.  condition1
             corresponds to one condition in the one way linear model,
             condition2 to the other.
"""

def initialize_database(db):
    """Set up the schema for SQLite3 handle *db*.

    """
    db.execute("""
               create table sample_group (
                   id integer primary key,
                   label text unique,
                   is_control boolean
               )""")
    db.execute("""
               create table samples (
                   id integer primary key,
                   sample_group integer references sample_group(id),
                   filename text,
                   n_reads integer
               )""")
    db.execute("""
               create table transcripts (
                  id integer primary key,
                  label text,
                  length integer
               )
               """)
    db.execute("""
               create table leftsites (
                   sample integer references samples(id),
                   transcript integer references transcripts(id),
                   position integer not null,
                   n integer not null default 0,
                   primary key (sample,transcript,position)
               )
               """)
    db.execute("""
               create table multiplicities (
                   id integer primary key,
                   sample integer references leftsites(sample),
                   n integer
               )
               """)
    db.execute("""
               create table multiplicity_entries (
                   id integer primary key,
                   transcript integer references leftsites(transcript),
                   position integer references leftsites(position),
                   multiplicity integer references multiplicities(id)
               )
               """)
    db.commit()
    vmsg("Initialized database.")

def insert_sample_group(db, label, is_control):
    db.execute("""insert into sample_group (label,is_control)
                  values (?,?)""", (label, is_control))
    (sample_group,) = db.execute("""select last_insert_rowid()""").fetchone()
    vmsg("Assigned sample group id %d to group with label %s" % (sample_group, label))
    return sample_group

def insert_sample(db, filename, sample_group):
    db.execute("""insert into samples (filename,sample_group)
                  values (?,?)""", (filename, sample_group))
    (sample,) = db.execute("""select last_insert_rowid()""").fetchone()
    vmsg("Assigned sample id %d" % sample)
    return sample

def insert_or_check_transcripts(db, sample, transcripts):
    if db.execute("""select count(id)>0 
                     from transcripts""").fetchone()[0] == 1:
        # Another call to load_sam has already loaded the transcripts
        # for this analysis.  Just check that the transcripts in this
        # file match those already loaded.
        vmsg("Transcripts already loaded.  Checking for integrity.")
        for i,h in enumerate(transcripts):
            q = db.execute("""select label,length from transcripts
                              where id=?""", (i,)).fetchone()
            if q == None:
                raise ValueError(("Failed checking transcripts against " + \
                                     "database: transcript in position " + \
                                     "%d with label %s does not exist in " + \
                                     "database.") % (i,h['SN']))
            else:
                (label,length) = q
                if label != h['SN'] or length != (h['LN']-38):
                    raise ValueError(("Transcript at position %d of %s does " + \
                                         "not match existing database.  " + \
                                         "Database had label %s with " + \
                                         "length %d; file had label %s " + \
                                         "with length %d.") % (i,filename,
                                                              h['SN'],h['LN'],
                                                              label,length+38))
    else:
        # The database has no transcripts.  Insert them.
        vmsg("Loading transcripts into database.")
        for i,h in enumerate(transcripts):
            db.execute("""insert into transcripts(id,label,length)
                          values (?,?,?)""", (i,h['SN'],h['LN']-38))

    vmsg("Initializing empty leftsites arrays.")
    for i,h in enumerate(transcripts):
        for p in range(h['LN']-38):
            db.execute("""insert into leftsites(sample,transcript,position,n) 
                          values (?,?,?,0)""", (sample,i,p))


def insert_reads_and_multiplicities(db, sample, samfile):
    n_reads = 0
    for readset in split_by_readname(samfile):
        n_reads += 1
        if len(readset) > 1:
            targets = tuple([(r.rname,r.pos) for r in readset])
            mid = (sample,targets).__hash__()
            if db.execute("""select id from multiplicities where id=?""", (mid,)).fetchone() == None:
                db.execute("""insert into multiplicities(id,sample,n)
                              values (?,?,1)""", (mid, sample))
                for (t,p) in targets:
                    db.execute("""insert into multiplicity_entries 
                                  (transcript,position,multiplicity) 
                                  values (?,?,?)""", (t,p,mid))
            else:
                db.execute("""update multiplicities set n=n+1
                              where id=?""", (mid,))
        for r in readset:
            db.execute("""update leftsites set n=n+1 where
                          sample=? and transcript=? and position=?""",
                       (sample,r.rname,r.pos))
    return n_reads


def load_sam(db, filename, sample_group):
    vmsg("Loading reads from %s with covariate %f" % (filename,sample_group))
    s = pysam.Samfile(filename)

    sample = insert_sample(db, filename, sample_group)
    insert_or_check_transcripts(db, sample, s.header['SQ'])
    n_reads = insert_reads_and_multiplicities(db, sample, s)
    db.execute("""update samples set n_reads=? where id=?""",
               (n_reads, sample))
            
    db.commit()
    s.close()
    return sample


def split_by_readname(samfile):
    """Return an iterator over the reads in *samfile* grouped by read name.

    The SAM file produced by bowtie is sorted by read name.  Often we
    want to work with all of the alignments of a particular read at
    once.  This function turns the flat list of reads into a list of
    lists of reads, where each sublist has the same read name.
    """
    last_read = None
    for r in samfile:
        if r.qname != last_read:
            if last_read != None:
                yield accum
            accum = [r]
            last_read = r.qname
        else:
            accum.append(r)
    yield accum    

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
            opts, args = getopt.getopt(argv, "hvl:", ["help"])
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
            if o in ("-l", ):
                try:
                    state.read_length = int(a)
                except ValueError, v:
                    raise Usage("Read length must be an integer, found %s" % a)
            else:
                raise Usage("Unhandled option: " + o)
        if len(args) != 3:
            raise Usage("samfiles_to_sqlite.py takes exactly three arguments.")
        condition1_files = args[0].split(',')
        for f in condition1_files:
            if not(os.path.exists(f)):
                raise Usage("Input file %s does not exist." % f)

        condition2_files = args[1].split(',')
        for f in condition2_files:
            if not(os.path.exists(f)):
                raise Usage("Input file %s does not exist." % f)

        db_filename = args[2]
        if os.path.exists(db_filename):
            raise Usage("Database file %s already exists" % db_filename)

        db = sqlite3.connect(db_filename)
        initialize_database(db)
        control_sample_group = insert_sample_group(db, "Controls", True)
        for f in condition1_files:
            load_sam(db, f, control_sample_group)
        other_sample_group = insert_sample_group(db, "Others", False)
        for f in condition2_files:
            load_sam(db, f, other_sample_group)


    
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

