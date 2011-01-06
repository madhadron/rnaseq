
import numpy as np
cimport numpy as np
from pymc import *
import sqlite3
from bag import *

def samples_of_group(db, sample_group):
    """Fetches the samples and numbers of reads in 'sample_group'

    Returns a dictionary with sample ID as keys and the number of
    reads in that sample as values.
    """
    c = db.execute("""select id,n_reads from samples
                      where sample_group = ?""", (sample_group,))
    samples = dict(c.fetchall())
    if samples == {}:
        raise ValueError("No samples associated to a group %d" % sample_group)
    else:
        return samples

def get_sample(db, sample_id, transcripts):
    """Fetch leftsites and multiplicities for 'transcripts' in 'sample_id'.

    'db' should be an SQLite3 database handle containing the sample
    groups.  'sample_id' is an integer telling which sample to fetch.
    Only those transcripts in 'transcripts' (a list of integers) will
    be returned.

    Returns a dictionary with the transcript IDs as keys, and a
    dictionary with key 'leftsites' referring to a numpy array and
    'multiplicities' referring to a Bag (see bag.py) of multiplicities
    as tuples of transcript IDs.
    """
    r = {}
    for t in transcripts:
        c = db.execute("""select n from leftsites where sample=? 
                          and transcript=? order by position asc""",
                       (sample_id, t))
        leftsites = np.array([x for (x,) in c])
        if len(leftsites) == 0:
            raise ValueError("No transcript with ID %d for sample %d in database" % (t,sample_id))
        multiplicities = Bag()
        c = db.execute("""select a.multiplicity, a.transcript, b.position 
                          from multiplicity_entries as a
                          join multiplicity_entries as b
                          on b.transcript = ? and 
                          a.transcript != ? and
                          a.multiplicity = b.multiplicity
                          join multiplicities as c
                          on a.multiplicity = c.id
                          where c.sample = ?""",
                       (t,t,sample_id))
        for m in group_by_first(c):
            pos = m[0][1]
            targets = tuple([x[0] for x in m])
            multiplicities.update([(pos,targets)])
        r[t] = {'leftsites': leftsites, 'multiplicities': multiplicities}
    return r

def alpha(float logp):
    """Calculates first parameter of Beta distribution.

    'logp' should be a double.
    """
    cdef double eta
    eta = np.exp(4.9 - 0.7*logp) - 1
    return np.exp(logp)*eta

def beta(float logp):
    """Calculates the second parameter of the Beta distribution.
    
    'log' can be a floating point number or a numpy array.
    """
    cdef double eta
    eta = np.exp(4.9 - 0.7*logp) - 1
    return (1 - np.exp(logp))*eta

def maintain_positive_parameters(minusmu=None, a=None):
    """Guard against feeding negative parameters to Beta distribution.

    Unconstrainted, 'minusmu' is Gamma distributed and 'a' is Cauchy
    distributed, so their combination can produce negative values as
    the parameters of the Beta distribution that gives the values of r
    below.  This function is used in a potential to prevent such
    cases.
    """
    if alpha(-1*minusmu + 0.5*a) <= 0 or \
            alpha(-1*minusmu - 0.5*a) <= 0 or \
            beta(-1*minusmu + 0.5*a) <= 0 or \
            beta(-1*minusmu - 0.5*a) <= 0:
        return -np.inf
    else:
        return 0


def multiplicity_correction(int T, int L, int transcript, dict rs, multiplicities):
    cdef np.ndarray[np.double_t, ndim=1] pm
    cdef double thisr, z
    cdef int position, multiplicity
    cdef tuple targets
    thisr = rs[transcript]
    pm = (thisr*T/L) * np.ones(L)
    for (position,targets),multiplicity in multiplicities.itercounts():
        z = sum(rs[k] for k in targets)
        pm[position] += multiplicity * z / (z + thisr)
    return pm

def make_observation(group_id, sample_id, transcript, rs, T, leftsites, multiplicities):
    """Creates multiread corrected random variables.

    Returns a list of random variables.  The first is a Deterministic
    which calculates a mean for the Poisson distribution with
    multiread corrections in 'multiplicities'.  The second is a
    Poisson distribution using that mean and associated with the
    observations in 'leftsites'.
    """
    L = len(leftsites)
    def _pm(rs = None):
        return multiplicity_correction(T, L, transcript, rs, multiplicities)
    pm_name = 'poisson_mean'+str(transcript)+'-group'+str(group_id)+'-'+str(sample_id)
    poisson_mean = Deterministic(eval=_pm,
                                 doc='Corrected mean for Poisson distribution',
                                 name=pm_name,
                                 parents={'rs': rs},
                                 trace=False, 
                                 verbose=0, 
                                 dtype=float,
                                 plot=False, 
                                 cache_depth=2)
    observation = Poisson('d'+str(transcript)+'-group'+str(group_id)+\
                              '-'+str(sample_id),
                          mu = poisson_mean,
                          observed=True, value=leftsites)
    return [poisson_mean, observation]


def build_model(db, group1, group2, transcripts):
    n_transcripts = db.execute("""select count(id) from transcripts""").fetchone()[0]
    n_reads = {1: samples_of_group(db, group1),
               2: samples_of_group(db, group2)}
    samples1 = n_reads[1].keys()
    samples2 = n_reads[2].keys()
    data = {1: dict([(s,get_sample(db, s, transcripts))
                     for s in n_reads[1].keys()]),
            2: dict([(s,get_sample(db, s, transcripts))
                     for s in n_reads[2].keys()])}
    [a,minusmu,maintain_beta,alphas,betas,r,d] = [{},{},{},{},{},{},{}]
    for t in transcripts:
        # t gets reassigned at each iteration, not redefined, so it
        # will be dynamically scoped if we use it as a free variable
        # in defining the model.  We assign tr the value of t, so it
        # will be properly lexically scoped, and use that instead.
        tr = t
        minusmu[t] = Gamma('minusmu'+str(t), 
                           alpha=5230.0/np.sqrt(n_transcripts),
                           beta=1/(2.1e-3 * np.sqrt(n_transcripts)))
        a[t] = Cauchy('a'+str(t), 0, 29)
        a[t].value = 0
        maintain_beta[t] = Potential(logp = maintain_positive_parameters,
                                     name = 'maintain_beta' + str(t),
                                     parents = {'minusmu': minusmu[t],
                                                'a': a[t]},
                                     doc = 'Maintain beta parameters positive',
                                     verbose = 0,
                                     cache_depth = 2)

    def make_alpha(cov):
        def _f(minusmu=None, a=None):
            return alpha(-1*minusmu + cov*a)
        return _f

    def make_beta(cov):
        def _f(minusmu=None, a=None):
            return beta(-1*minusmu + cov*a)
        return _f


    # All group 1 samples use a covariate of 0.5, all group 2
    # samples a covariate of -0.5.
    r[1] = {}
    alphas[1] = {}
    betas[1] = {}
    for sample in samples1:
        s = sample
        r[1][s] = {}
        alphas[1][s] = {}
        betas[1][s] = {}
        for t in transcripts:
            tr = t
            alphas[1][s][t] = Deterministic(eval=make_alpha(0.5),
                                            doc='',
                                            name='alpha'+str(t)+'-group1-'+str(s),
                                            parents={'minusmu':minusmu[tr],
                                                     'a':a[tr]},
                                            trace=False,
                                            verbose=0,
                                            dtype=float,
                                            plot=False,
                                            cache_depth=2)
            betas[1][s][t] = Deterministic(eval=make_beta(0.5),
                                          doc='',
                                          name='beta'+str(t)+'-group1-'+str(s),
                                          parents={'minusmu':minusmu[tr],
                                                   'a':a[tr]},
                                          trace=False,
                                          verbose=0,
                                          dtype=float,
                                          plot=False,
                                          cache_depth=2)
            r[1][s][t] = Beta('r'+str(t)+'-group1-'+str(s),
                              alpha=alphas[1][s][tr],
                              beta=betas[1][s][tr],
                              trace = True)
    r[2] = {}
    alphas[2] = {}
    betas[2] = {}
    for sample in samples2:
        s = sample
        r[2][s] = {}
        alphas[2][s] = {}
        betas[2][s] = {}
        for t in transcripts:
            tr = t
            alphas[2][s][t] = Deterministic(eval=make_alpha(-0.5),
                                            doc='',
                                            name='alpha'+str(t)+'-group2-'+str(s),
                                            parents={'minusmu':minusmu[tr],
                                                     'a':a[tr]},
                                            trace=False,
                                            verbose=0,
                                            dtype=float,
                                            plot=False,
                                            cache_depth=2)
            betas[2][s][t] = Deterministic(eval=make_beta(-0.5),
                                          doc='',
                                          name='beta'+str(t)+'-group2-'+str(s),
                                          parents={'minusmu':minusmu[tr],
                                                   'a':a[tr]},
                                          trace=False,
                                          verbose=0,
                                          dtype=float,
                                          plot=False,
                                          cache_depth=2)
            r[2][s][t] = Beta('r'+str(t)+'-group2-'+str(s),
                              alpha=alphas[2][s][tr],
                              beta=betas[2][s][tr],
                              trace = True)

    d[1] = {}
    for sample in samples1:
        s = sample
        d[1][s] = []
        for t in transcripts:
            d[1][s].extend(make_observation(1, s, t, r[1][s], n_reads[1][s], 
                                            data[1][s][t]['leftsites'],
                                            data[1][s][t]['multiplicities']))
    d[2] = {}
    for sample in samples2:
        s = sample
        d[2][s] = []
        for t in transcripts:
            d[2][s].extend(make_observation(2, s, t, r[2][s], n_reads[2][s], 
                                            data[2][s][t]['leftsites'],
                                            data[2][s][t]['multiplicities']))
            
    return MCMC([minusmu, a, maintain_beta, alphas, betas, r, d])



