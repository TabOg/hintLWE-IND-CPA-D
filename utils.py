from estimator import *
from functools import partial
from sage.all import oo, log

# this function returns the bit security according to the cost model, reduced basis shape, and attacks considered in the HE Standard.
# See https://github.com/gong-cr/FHE-Security-Guidelines/ for further details.
def HE_standard_LWE_hardness(logn, logQ, secret_sigma, error_sigma=3.19):
    attacks = [
        partial(LWE.primal_usvp, red_cost_model=RC.MATZOV),
        partial(LWE.dual_hybrid, red_cost_model=RC.MATZOV),
        partial(LWE.primal_bdd, red_cost_model=RC.MATZOV),
    ]
    
    # when the parameters are small, we also consider the primal hybrid with the following settings
    primal_hybrid = partial(LWE.primal_hybrid, mitm=False, babai=False, red_cost_model=RC.MATZOV)
    
    params = LWE.Parameters(2 ** logn, 2 ** logQ, Xs=ND.DiscreteGaussian(secret_sigma), Xe=ND.DiscreteGaussian(error_sigma), m=oo)
    
    # we loop through all the attacks, and return the lowest estimate in bits
    min_cost = oo
    for attack in attacks:
        try:
            cost = attack(params=params)
            if log(cost["rop"], 2).n() < min_cost:
                min_cost = log(cost["rop"], 2).n()
        except:
            continue
    
    if logn <= 14:
        try:
            cost = primal_hybrid(params=params)
            if log(cost["rop"], 2).n() < min_cost:
                min_cost = log(cost["rop"], 2).n()
        except:
            pass
    return min_cost