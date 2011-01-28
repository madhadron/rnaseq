"""
workflow.py

Python script which hooks into the HTSStation frontend to do analysis
of RNASeq data.

"""
import getopt
import os
import sys
import pysam
import sqlite3
from bbcflib import *
from bein.util import *
from rnaseq import *

usage = """workflow.py [-vh] [-l readlen] working_lims config_lims job_key

-v           Run verbosely
-h           Print this message and exit
-l readlen   Reads have length 'readlen' in SAM/BAM files
working_lims MiniLIMS where RNASeq executions and files will be stored.
config_lims  MiniLIMS containing a pickled ConfigParser under the alias 'config'
job_key      Alphanumeric key specifying the job
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class State(object):
    def __init__(self):
        self.verbose = False
        self.read_length = 38
        self.working_lims = None
        self.config_lims = None
        self.job_key = None
        self.config = None

state = State()

def vmsg(msg):
    if state.verbose:
        print >>sys.stderr, msg

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        try:
            opts, args = getopt.getopt(argv, "hvl:", 
                                       ["help","read-length"])
        except getopt.error, message:
            raise Usage(message)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-v", ):
                state.verbose=True
                vmsg("Running verbosely.")
            elif o in ("-l", "--read-length"):
                try:
                    state.read_length = int(a)
                    vmsg("Using read length %d" % state.read_length)
                except ValueError, v:
                    raise Usage("Read length must be an integer, found %s" % a)
            else:
                raise Usage("Unhandled option: " + o)
        if len(args) != 3:
            raise Usage("samfiles_to_sqlite.py takes exactly three arguments.")
        try:
            state.job_key = args[2]
            vmsg("Job key is %s" % args[2])
        except ValueError, v:
            raise Usage("job_id must be an integer, found %s", str(args[2]))

        state.working_lims = MiniLIMS(args[0])
        vmsg("Connected to %s as working LIMS" % args[0])

        if not(os.path.exists(args[1])):
            raise Usage("config_lims %s does not exist." % args[1])
        else:
            state.config_lims = MiniLIMS(args[1])
            state.config = use_pickle(state.config_lims, "config")
            vmsg("Loaded configuration from LIMS %s" % args[1])

        # Fetch job information from frontend
        frontend = Frontend('http://htsstation.vital-it.ch/rnaseq/')
        try:
            job = frontend.job(state.job_key)
            vmsg("Fetched job from frontend at %s" % 'http://htsstation.vital-it.ch/rnaseq/')
        except TypeError, t:
            raise Usage("No such job with key %s at frontend %s" % (state.job_key,
                            'http://htsstation.vital-it.ch/rnaseq/'))                                        
        # Get bowtie index path from GenRep
        genrep = GenRep('http://bbcftools.vital-it.ch/genrep/',
                        '/scratch/frt/yearly/genrep/nr_assemblies/cdna_bowtie')
        assembly = genrep.assembly(job.assembly_id)
        vmsg("Fetched assembly %d from GenRep" % job.assembly_id)

        # Build database for inference from DAF LIMS files.
        # Fetch all the FASTQ files from the DAF LIMS, run bowtie on
        # them, then load_sam on each of them.  Bowtie is run in
        # parallel, but load_sam unfortunately has to be done in
        # sequence due to database locking issues.
        daflims = DAFLIMS(username='jrougemont', password='cREThu6u')
        vmsg("Connected to DAF LIMS")
        with execution(state.working_lims) as ex:
            db_name = unique_filename_in()
            db = sqlite3.connect(db)
            initialize_database(db)
            fastqfiles = {}
            for gid,g in job.groups.iteritems():
                fastqfiles[gid] = {}
                insert_sample_group(db, g.name, g.control, gid)
                for rid,r in g['runs'].iteritems():
                    filename = unique_filename_in()
                    fastqfiles[gid][rid] = daflims.fetch(r.facility, r.machine, 
                                                         r.run, r.lane, filename)
                    vmsg("Fetched %s/%s/%d/%d from DAF LIMS as %s" %
                         (r.facility, r.machine, r.run, r.lane, filename))
            def _align(filename):
                return background(parallel_bowtie_lsf, 
                                  ex, assembly.index_path,
                                  filename, bowtie_args="-Sqa",
                                  add_nh_flags=True)
            samfile_futures = deepmap(_align, fastqfiles)
            samfiles = deepmap(lambda q: q.wait(), samfile_futures)
            for gid in fastqfiles.keys():
                [load_sam(db, f, gid) for f in fastqfiles[gid].itervalues()]
            ex.add(db_name)
            vmsg("Finished loading SAM files into database.")

        # Find separable subproblems
        subproblems = find_subproblems(db)

        # Make list of pairs of groups
        # If all groups are control or all are not control, do every
        # combination of two groups.  Otherwise do every combination
        # of a control and a non-control group.
        if all([x['control'] for x in job.groups.itervalues()]) or \
                not(any([x['control'] for x in job.groups.itervalues()])):
            # case: all are control or none are control
            pairs = [(x,y) for x in job.groups.iterkeys()
                     for y in job.groups.iterkeys()
                     where x != y]
            vmsg("Going to run all against all: %s" % str(pairs))
        else:
            # case: mixed
            control_group_ids = [gid for gid,g in job.groups.iteritems()
                                 if g['control'] == True]
            other_group_ids = [gid for gid,g in job.groups.iteritems()
                               if g['control'] == False]
            pairs = [(x,y) for x in control_group_ids for y in other_group_ids]
            vmsg("Going to run control against non-control: %s" % str(pairs))

        # Run parallel executions which run each subproblem in parallel on a pair
        # leaving a pickle file, then run another script to add all the pickles to
        # to the database, and write a summary file which it adds to the MiniLIMS.

        # Send a report email of the run

        # Call back to the frontend that the job is complete

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())

