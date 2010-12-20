import numpy as np
cimport numpy as np
from pymc import *

def build_model(db, transcripts):
    """Leftsite arrays is now a list of lists of arrays."""
    covariates = [x for (x,) in db.execute("""select distinct covariate 
                                              from samples""")]
    n_transcripts = db.execute("""select count(id) from transcripts""").fetchone()[0]

    def maintain_positive_parameters(minusmu=None, a=None):
        def logp(c):
            return -1*minusmu + a*c
        def eta(c):
            return np.exp(4.9 - 0.7*logp(c)) - 1
        def p(c):
            return np.exp(logp(c))
        if any([p(c) <= 0 for c in covariates]) or \
                any([1-p(c) <= 0 for c in covariates]):
            return -inf
        else:
            return 0

    def make_logp(c):
        def _logp(minusmu=None, a=None):
            return -1*minusmu + a*c
        return _logp

    def make_eta(c):
        def _eta(logp=None):
            return np.exp(4.9 - 0.7*logp) - 1
        return _eta

    def make_p(c):
        def _p(logp=None):
            return np.exp(logp)
        return _p

    [a,minusmu,maintain_beta,logp,eta,p,r,poisson_mean] = [{},{},{},{},{},{},{},{}]
    for t in transcripts:
        minusmu[t] = Gamma('minusmu'+str(t), alpha=11, beta=43.5/np.sqrt(n_transcripts))
        a[t] = Cauchy('a'+str(t), 0, 29)
        a[t].value = 0
        maintain_beta[t] = Potential(logp = maintain_positive_parameters,
                                     name = 'maintain_beta' + str(t),
                                     parents = {'minusmu': minusmu[t],
                                                'a': a[t]},
                                     doc = 'Maintain beta parameters positive',
                                     verbose = 0,
                                     cache_depth = 2)

        for c in covariates:
            # c gets reassigned at each iteration, not redefined, so
            # it will by dynamically scoped if used as a free variable
            # in lambda.  By assigning cv, we get a lexically scoped
            # free variable.
            cv = c 
            logp[(t,c)] = Deterministic(eval=make_logp(cv),
                                        name='logp'+str(t)+str(cv),
                                        parents={'minusmu': minusmu[t],
                                                 'a': a[t]},
                                        doc='log of p',
                                        trace=True,
                                        verbose=0,
                                        dtype=float,
                                        plot=False,
                                        cache_depth=2)
            eta[(t,c)] = Deterministic(eval=make_eta(cv),
                                       name='eta'+str(t)+str(cv),
                                       parents={'logp': logp[(t,cv)]},
                                       doc='eta',
                                       trace=True,
                                       verbose=0,
                                       dtype=float,
                                       plot=False,
                                       cache_depth=2)
            p[(t,c)] = Deterministic(eval=make_p(cv),
                                     name='p'+str(t)+str(cv),
                                     parents={'logp': logp[(t,cv)]},
                                     doc='p',
                                     trace=True,
                                     verbose=0,
                                     dtype=float,
                                     plot=False,
                                     cache_depth=2)
            r[(t,c)] = Beta('r'+str(t)+str(cv),
                            alpha=p[(t,cv)]*eta[(t,cv)],
                            beta=(1-p[(t,cv)])*eta[(t,cv)])



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


    return MCMC([minusmu, a, maintain_beta, logp, eta, p, r])
