### This file corresponds to section 6.1

## all outputs are formatted with Figure 1 in mind

import math
from utils import *
from sage.all import pi, sqrt, floor

# this function returns the bit security of a set of parameters
def original_bit_security_levels(parameters):
    security_levels = []
    for (logn, logq, sigma_s, sigma_e) in parameters:
        security_level = HE_standard_LWE_hardness(logn, logq, sigma_s, sigma_e)
        security_levels.append((logn,security_level))
    return security_levels

# the security of HintLWE (Decision) via the reduction from LWE of Corollary 1.
# we report the security of the tightest possible reduction
def hint_lwe_decision_security_levels(parameters, original_security):
    # store the derived security level
    kappa_primes = []
    # store the normalised flooding standard deviation in bits
    normalised_flooding_stddevs = []
    
    # this is the variance of the secret after hints, from Corollary 1.
    # Note that we have a circular dependence between the security level and the secret variance
    sigma_prime = lambda logn, kappa_prime: (log(4) + logn * log(2) + log(2 ** ((kappa_prime + 2) / 2) + 1)) / pi ** 2
    
    for (logn, logq, sigma_s, sigma_e), (_, original_kappa) in zip(parameters, original_security):
        found_reduction = False
        # we're looking for a kappa_prime such that sigma_prime(kappa_prime) < sigma_s ** 2 / 2 
        # AND LWE with secret of width sigma_prime is kappa bit secure, with kappa_prime <= kappa - 4
        for kappa_prime in reversed(range(floor(original_kappa))):
            sigma_prime_ = sigma_prime(logn, kappa_prime)
            # secret is not wide enough to accomodate reduction to this kappa_prime
            foo = sigma_s ** 2 / 2
            if sigma_prime_.n() >= sigma_s ** 2 / 2:
                continue
            # now we check if this LWE instance is kappa bit secure
            kappa = HE_standard_LWE_hardness(logn, logq, sigma_prime_, sigma_e)
            if kappa_prime > kappa - 4:
                continue
            
            kappa_primes.append((logn, kappa_prime))
            flooding_var = 1 / (1 / (2 * sigma_prime_) - 1 / sigma_s ** 2)
            normalised_flooding_stddevs.append((logn, 0.5 * log(flooding_var, 2).n()))
            found_reduction = True
            break
        if not found_reduction:
            print(f"parameters {(logn, logq, sigma_s, sigma_e)} do not accomodate reduction from LWE: try increasing sigma_s")
    return kappa_primes, normalised_flooding_stddevs
    
# the security of HintLWE (Search) via the reduction from LWE of Corollary 2. This corresponds to the red line in Fig 1.b.
# we report the security of the tightest possible reduction, by searching over 1 to l_max
# to pipe into later results, we also report sigma' squared
def hint_lwe_search_security_levels(parameters, l_max=128):
    # store the derived security level
    kappa_primes = []
    # store the normalised flooding standard deviation in bits
    normalised_flooding_stddevs = []
    
    # this is the variance of the secret after hints, from Corollary 2
    sigma_prime = lambda logn, l: (log(4) + logn * log(2) - log(1 - 2 ** (-l / 2))) / pi ** 2
    for (logn, logq, sigma_s, sigma_e) in parameters:
        found_reduction = False
        # tightest possible reduction: l = 1
        for l in range(1, l_max + 1):
            sigma_prime_ = sigma_prime(logn, l).n()
            
            # check if the original secret is wide enough to accomodate reduction
            if sigma_prime_ >= sigma_s ** 2 / 2:
                continue
            kappa_prime = HE_standard_LWE_hardness(logn, logq, sqrt(sigma_prime_), sigma_e) - l
            kappa_primes.append((logn, kappa_prime))
            flooding_var = 1 / (1 / (2 * sigma_prime_) - 1 / sigma_s ** 2)
            normalised_flooding_stddevs.append((logn, 0.5 * log(flooding_var, 2).n()))
            found_reduction = True
            break
        if not found_reduction:
            print(f"parameters {(logn, logq, sigma_s, sigma_e)} do not accomodate reduction from LWE: increase sigma_s or l_max")
    return kappa_primes, normalised_flooding_stddevs

# how much noise flooding is required for the decision reduction to give that the HintLWE problem with the given parameters is `security level` bit secure.
def hint_lwe_decision_normalised_noise_flooding(parameters, security_level):
    normalised_flooding_stddevs = []
    sigma_prime = lambda logn, security_level: (log(4) + logn * log(2) + log(2 ** ((security_level + 2) / 2) + 1)) / pi ** 2
    for (logn, logq, sigma_s, sigma_e) in parameters:
        sigma_prime_ = sigma_prime(logn, security_level)
        if sigma_prime_ > 0.5 * sigma_s ** 2:
            # can't achieve this security level with these parameter levels: continue
            continue
        else:
            # if LWE with secret stddev sigma_prime is kappa bits secure, and security_level <= kappa - 4, done
            kappa = HE_standard_LWE_hardness(logn, logq, sqrt(sigma_prime_), sigma_e)
            if security_level <= kappa - 4:
                flooding_var = 1 / (1 / (2 * sigma_prime_) - 1 / sigma_s ** 2)
                normalised_flooding_stddevs.append((logn, 0.5 * log(flooding_var, 2).n()))
    return normalised_flooding_stddevs
        
def hint_lwe_search_normalised_noise_flooding(parameters, security_level):
    normalised_flooding_stddevs = []
    sigma_prime = lambda logn, l: (log(4) + logn * log(2) - log(1 - 2 ** (-l / 2))) / pi ** 2
    for (logn, logq, sigma_s, sigma_e) in parameters:
        highest_sigma_prime = sigma_prime(logn, 1)
        highest_kappa = HE_standard_LWE_hardness(logn, logq, sqrt(highest_sigma_prime), sigma_e)
        if (security_level + 1 - highest_kappa) > -10 ** (-10):
            # can't achieve this security level with these parameter levels: continue
            continue
        
        lowest_sigma_prime = (log(4) + logn * log(2)) / pi ** 2
        if lowest_sigma_prime > 0.5 * sigma_s ** 2:
            # can't achieve this security level with these parameter levels: continue
            continue
        lowest_kappa = HE_standard_LWE_hardness(logn, logq, sqrt(lowest_sigma_prime), sigma_e)
        for l in reversed(range(1, int(math.ceil(lowest_kappa - security_level)) + 1)):
            sigma_prime_ = sigma_prime(logn, l)
            if sigma_prime_ > 0.5 * sigma_s ** 2:
                continue
            kappa = HE_standard_LWE_hardness(logn, logq, sqrt(sigma_prime_), sigma_e)
            kappa_prime = kappa - l
            if (kappa_prime - security_level) >= -10 ** (-10):
                flooding_var = 1 / (1 / (2 * sigma_prime_) - 1 / sigma_s ** 2)
                normalised_flooding_stddevs.append((logn, 0.5 * log(flooding_var, 2).n()))
                break
    return normalised_flooding_stddevs
    
if __name__ == "__main__":
    # this script generates the data for Fig 1 and 2
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
    
    original_security_levels=[(10, 131.903812370381), (11, 128.744594166727), (12, 128.934036200343), (13, 128.193639623478), (14, 128.334392294043), (15, 128.078948557141), (16, 128.021734664213), (17, 128.126809683405)]
    
    decision_kappa_primes, decision_normalised_flooding_stddevs = hint_lwe_decision_security_levels(parameters, original_security_levels)
    # the red line in Fig 1.a.
    print(f"{decision_kappa_primes=}")
    print(f"{decision_normalised_flooding_stddevs=}")
    print()
    
    decision_normalised_flooding_stddev = hint_lwe_decision_normalised_noise_flooding(parameters, 118)
    print(decision_normalised_flooding_stddev)
    search_kappa_primes, search_normalised_flooding_stddevs = hint_lwe_search_security_levels(parameters)
    # the red line in Fig 1.b.
    print(f"{search_kappa_primes=}")
    print(f"{search_normalised_flooding_stddevs=}")
    print()
    
    # # data for Fig 2.a. and 2.b.
    for security_level in [120, 100, 80]:
        decision_normalised_noise_flooding_stddev = hint_lwe_decision_normalised_noise_flooding(parameters, security_level)
        print(f"{security_level=}:\t{decision_normalised_noise_flooding_stddev=}")
    print()
    
    for security_level in [120, 100, 80]:
        search_normalised_noise_flooding_stddev = hint_lwe_search_normalised_noise_flooding(parameters, security_level)
        print(f"{security_level=}:\t{search_normalised_noise_flooding_stddev=}")
    
    
## output:
# original_security_levels=[(10, 131.903812370381), (11, 128.744594166727), (12, 128.934036200343), (13, 128.193639623478), (14, 128.334392294043), (15, 128.078948557141), (16, 128.021734664213), (17, 128.126809683405)]

# decision_kappa_primes=[(10, 118), (11, 116), (12, 114), (13, 112), (14, 110), (15, 108), (16, 106), (17, 104)]
# decision_normalised_flooding_stddevs=[(10, 5.33793402193495), (11, 5.33793402193497), (12, 5.33793402193495), (13, 5.33793402193495), (14, 5.33793402193495), (15, 5.33793402193497), (16, 5.33793402193495), (17, 5.33793402193495)]

# search_kappa_primes=[(10, 122.715786881075), (11, 123.669580650461), (12, 125.837449190988), (13, 126.198049311874), (14, 126.925504468014), (15, 126.808342322618), (16, 126.981094074453), (17, 126.886665052395)]
# search_normalised_flooding_stddevs=[(10, 0.628012845817682), (11, 0.690977607639125), (12, 0.750845729343373), (13, 0.808032621561283), (14, 0.862883213741055), (15, 0.915687528719480), (16, 0.966692221128444), (17, 1.01610927748317)]

# security_level=120:     decision_normalised_noise_flooding_stddev=[]
# security_level=100:     decision_normalised_noise_flooding_stddev=[(10, 3.54720653853858), (11, 3.80683151509908), (12, 4.19411064008535), (13, 5.02076343018411)]
# security_level=80:      decision_normalised_noise_flooding_stddev=[(10, 2.64093356399580), (11, 2.73167067411133), (12, 2.83016417567041), (13, 2.93851227203193), (14, 3.05971885919284), (15, 3.19831311367933), (16, 3.36160890148374), (17, 3.56260895206806)]

# security_level=120:     search_normalised_noise_flooding_stddev=[(10, 0.551653685990066), (11, 0.604723116640469), (12, 0.654934374219918), (13, 0.712981272652456), (14, 0.771835026098540), (15, 0.828139617447590), (16, 0.882218547085110), (17, 0.934344977829440)]
# security_level=100:     search_normalised_noise_flooding_stddev=[(10, 0.507263701259074), (11, 0.577023429178899), (12, 0.642705792317042), (13, 0.704918455758421), (14, 0.764139939850901), (15, 0.820764608820826), (16, 0.875123638412292), (17, 0.927496239692706)]
# security_level=80:      search_normalised_noise_flooding_stddev=[(10, 0.507227843676999), (11, 0.576999631793473), (12, 0.642694555115769), (13, 0.704910910612041), (14, 0.764132740962542), (15, 0.820757711200289), (16, 0.875117004346347), (17, 0.927489837197405)]

