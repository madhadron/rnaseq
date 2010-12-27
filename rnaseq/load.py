import sqlite3
import pysam

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
    db.execute("""
               create table inferences (
                   id integer primary key,
                   group1 integer references sample_group(id),
                   group2 integer references sample_group(id),
                   unique (group1,group2),
                   check (group1 <= group2)
               )
               """)
    db.execute("""
               create table posterior_samples (
                   inference integer references inferences(id),
                   transcript integer references transcripts(id),
                   variable text not null,
                   sample integer not null,
                   value float not null,
                   primary key (inference,transcript,variable,sample)
               )
               """)
    db.commit()

def insert_sample_group(db, label, is_control, group_id=None):
    if group_id != None:
        x = db.execute("""select id from sample_group where id=?""", (group_id,)).fetchall()
        if x != None:
            raise ValueError("Group %d already exists in database." % group_id)
    elif group_id == None:
        db.execute("""insert into sample_group (label,is_control)
                      values (?,?)""", (label, is_control))
    else:
        db.execute("""insert into sample_group (id,label,is_control)
                      values (?,?,?)""", (group_id, label, is_control))
    (sample_group,) = db.execute("""select last_insert_rowid()""").fetchone()
    return sample_group

def insert_sample(db, filename, sample_group):
    db.execute("""insert into samples (filename,sample_group)
                  values (?,?)""", (filename, sample_group))
    (sample,) = db.execute("""select last_insert_rowid()""").fetchone()
    return sample

def insert_or_check_transcripts(db, sample, transcripts):
    if db.execute("""select count(id)>0 
                     from transcripts""").fetchone()[0] == 1:
        # Another call to load_sam has already loaded the transcripts
        # for this analysis.  Just check that the transcripts in this
        # file match those already loaded.
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
        for i,h in enumerate(transcripts):
            db.execute("""insert into transcripts(id,label,length)
                          values (?,?,?)""", (i,h['SN'],h['LN']-38))

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
