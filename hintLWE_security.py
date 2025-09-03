### This file corresponds to section 6.1

## all outputs are formatted with Figure 1 in mind

from utils import *
from sage.all import pi, sqrt, floor

# this function returns the bit security of a set of parameters
def original_bit_security_levels(parameters):
    security_levels = []
    for (logn, logQ, sigma_s, sigma_e) in parameters:
        security_level = HE_standard_LWE_hardness(logn, logQ, sigma_s, sigma_e)
        security_levels.append((logn,security_level))
    return security_levels

# the security of HintLWE (Decision) via the reduction from LWE of Corollary 1.
# we report the security of the tightest possible reduction
def hint_lwe_decision_security_levels(parameters, original_security, l = 4):
    # store the derived security level
    kappa_primes = []
    # store the width of the secret after hints: needed to calculate flooding variance
    sigma_primes = []
    
    # this is the variance of the secret after hints, from Corollary 1.
    # Note that we have a circular dependence between the security level and the secret variance
    sigma_prime = lambda logn, kappa, l: (log(4) + logn * log(2) + log(2 ** ((kappa - l + 2) / 2) + 1)) / pi ** 2
    
    for (logn, logQ, sigma_s, sigma_e), (_, original_kappa) in zip(parameters, original_security):
        found_reduction = False
        # we're looking for a (kappa, l) such that sigma_prime(kappa, l) < sigma_s ** 2 / 2 AND LWE with secret of width sigma_prime is kappa bit secure
        # to make this run in a reasonable time, we fix l = 4 (the most efficient reduction) and try decreasing kappa
        for kappa in reversed(range(floor(original_kappa))):
            sigma_prime_ = sigma_prime(logn, kappa, l)
            
            # secret is not wide enough to accomodate reduction
            if sigma_prime_.n() >= sigma_s ** 2 / 2:
                continue
            
            # now we check if this LWE instance is kappa bit secure
            security = HE_standard_LWE_hardness(logn, logQ, sigma_prime_, sigma_e)
            if security < kappa:
                continue
            
            kappa_primes.append((logn, kappa - l))
            sigma_primes.append((logn, sigma_prime_.n()))
            found_reduction = True
            break
        if not found_reduction:
            print(f"parameters {(logn, logQ, sigma_s, sigma_e)} do not accomodate reduction from LWE: try increasing sigma_s")
    return kappa_primes, sigma_primes
    
# the security of HintLWE (Search) via the reduction from LWE of Corollary 2. This corresponds to the red line in Fig 1.b.
# we report the security of the tightest possible reduction, by searching over 1 to l_max
# to pipe into later results, we also report sigma' squared
def hint_lwe_search_security_levels(parameters, l_max=128):
    # store the derived security level
    kappa_primes = []
    # store the width of the secret after hints: needed to calculate flooding variance
    sigma_primes = []
    
    # this is the variance of the secret after hints, from Corollary 2
    sigma_prime = lambda logn, l: (log(4) + logn * log(2) - log(1 - 2 ** (-l / 2))) / pi ** 2
    for (logn, logQ, sigma_s, sigma_e) in parameters:
        found_reduction = False
        # tightest possible reduction: l = 1
        for l in range(1, l_max + 1):
            sigma_prime_ = sigma_prime(logn, l).n()
            
            # check if the original secret is wide enough to accomodate reduction
            if sigma_prime_ >= sigma_s ** 2 / 2:
                continue
            security_level = HE_standard_LWE_hardness(logn, logQ, sqrt(sigma_prime_), sigma_e) - l
            kappa_primes.append((logn, security_level))
            sigma_primes.append((logn, sigma_prime_))
            found_reduction = True
            break
        if not found_reduction:
            print(f"parameters {(logn, logQ, sigma_s, sigma_e)} do not accomodate reduction from LWE: increase sigma_s or l_max")
    return kappa_primes, sigma_primes

if __name__ == "__main__":
    # this script generates the data for Fig 1.a. and Fig 1.b. of the paper
    # in addition, we generate the corresponding sigma' squared values to calculate the flooding noise variance
    # these are the 128 bit secure parameters from Table 5.2 of the HE Community Guidelines, column "Gaussian".
    # (logn, logq, sigma_s, sigma_e)
    parameters = [
        (10, 28, 3.19, 3.19),
        (11, 55, 3.19, 3.19),
        (12, 108, 3.19, 3.19),
        (13, 216, 3.19, 3.19),
        (14, 432, 3.19, 3.19),
        (15, 870, 3.19, 3.19),
        (16, 1749, 3.19, 3.19),
        (17, 3525, 3.19, 3.19)
    ]
    
    # the blue lines in Fig 1.a. and 1.b.
    original_security_levels = original_bit_security_levels(parameters=parameters)
    print(f"# {original_security_levels=}")
    print()
    
    
    decision_kappa_primes, decision_sigma_primes = hint_lwe_decision_security_levels(parameters, original_security_levels)
    # the red line in Fig 1.a.
    print(f"{decision_kappa_primes=}")
    print(f"{decision_sigma_primes=}")
    print()
    
    search_kappa_primes, search_sigma_primes = hint_lwe_search_security_levels(parameters)
    # the red line in Fig 1.b.
    print(f"{search_kappa_primes=}")
    print(f"{search_sigma_primes=}")
    
    
## output:
# original_security_levels=[(10, 131.903812370381), (11, 128.744594166727), (12, 128.934036200343), (13, 128.193639623478), (14, 128.334392294043), (15, 128.078948557141), (16, 128.021734664213), (17, 128.126809683405)]

# decision_kappa_primes=[(10, 118), (11, 116), (12, 114), (13, 112), (14, 110), (15, 108), (16, 106), (17, 104)]
# decision_sigma_primes=[(10, 5.05659547963317), (11, 5.05659547963317), (12, 5.05659547963317), (13, 5.05659547963317), (14, 5.05659547963317), (15, 5.05659547963317), (16, 5.05659547963317), (17, 5.05659547963317)]

# search_kappa_primes=[(10, 122.715786881075), (11, 123.669580650461), (12, 125.837449190988), (13, 126.198049311874), (14, 126.925504468014), (15, 126.808342322618), (16, 126.981094074453), (17, 126.886665052395)]
# search_sigma_primes=[(10, 0.967182974726449), (11, 1.03741346749913), (12, 1.10764396027181), (13, 1.17787445304450), (14, 1.24810494581718), (15, 1.31833543858986), (16, 1.38856593136255), (17, 1.45879642413523)]


    

