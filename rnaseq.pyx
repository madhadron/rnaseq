import numpy as np
cimport numpy as np
from pymc import *

cpdef double alpha(double logp):
    """Calculates first parameter of Beta distribution.
    
    'logp' should be a double.
    """
    cdef double eta
    eta = np.exp(4.9 - 0.7*logp) - 1
    return np.exp(logp)*eta

cpdef double beta(double logp):
    """Calculates the second parameter of the Beta distribution.
    
    'log' can be a floating point number or a numpy array.
    """
    cdef double eta
    eta = np.exp(4.9 - 0.7*logp) - 1
    return (1 - np.exp(logp))*eta

cpdef double maintain_positive_parameters(double minusmu=None, double a=None):
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
        return -inf
    else:
        return 0

def build_model(db, group1, group2):
    n_transcripts = db.execute("""select count(id) from transcripts""").fetchone()[0]


    [a,minusmu,maintain_beta,r] = [{},{},{},{}]
    for t in transcripts:
        # t gets reassigned at each iteration, not redefined, so it
        # will be dynamically scoped if we use it as a free variable
        # in defining the model.  We assign tr the value of t, so it
        # will be properly lexically scoped, and use that instead.
        tr = t
        minusmu[t] = Gamma('minusmu'+str(t), alpha=11, beta=43.5/np.sqrt(n_transcripts),
                           trace=True)
        a[t] = Cauchy('a'+str(t), 0, 29, trace=True)
        a[t].value = 0
        maintain_beta[t] = Potential(logp = maintain_positive_parameters,
                                     name = 'maintain_beta' + str(t),
                                     parents = {'minusmu': minusmu[t],
                                                'a': a[t]},
                                     doc = 'Maintain beta parameters positive',
                                     verbose = 0,
                                     cache_depth = 2,
                                     trace = False)

        # All group 1 samples use a covariate of 0.5, all group 2
        # samples a covariate of -0.5.
        for sample in group1:
            s = sample
            r[(t,1,s)] = Beta('r'+str(t)+'-group1-'+s,
                            alpha=alpha(-1*minusmu[tr] + 0.5*a[tr]),
                            beta=beta(-1*minusmu[tr] + 0.5*a[tr]),
                            trace = False)

        for sample in group2:
            r[(t,2,s) = Beta('r'+str(t)+'-group2-'+s,
                             alpha=alpha(-1*minusmu[tr] - 0.5*a[tr]),
                             beta=beta(-1*minusmu[tr] - 0.5*a[tr]),
                             trace = False)

            # poisson_mean[(t,c)] = Deterministic(eval=make_poisson_mean(i),
            #                                     name='poisson_mean'+str(t)+str(c),
            #                                     parents={'r': r},
            #                                     doc='Corrected Poisson mean',
            #                                     trace = True,
            #                                     verbose = 0,
            #                                     dtype = float,
            #                                     plot = False,
            #                                     cache_depth = 2)

    # ds = [Poisson('d'+names[i], mu=poisson_means[i],
    #               observed=True, value=mapping[i]['leftsites'])
    #       for i in range(N)]


    return MCMC([minusmu, a, maintain_beta, r])
