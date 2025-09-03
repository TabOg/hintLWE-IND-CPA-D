# this code corresponds to section 6.2
from sage.all import log, pi, sqrt

# this is the required flooding noise from the prior art as a function of a total noise bound E, number of decryptions, and a target security level
def prior_flooding_noise(t, logE, target_security):
    total_flooding_noise = logE + target_security / 2 + log(64 * t, 2) / 2
    return total_flooding_noise.n()

# this is the required flooding noise from our Theorem 3
def bit_security_flooding_noise(logE, target_security):
    total_flooding_noise = logE + target_security / 2 + log(4, 2) / 2
    return total_flooding_noise.n()

# our flooding noise is linear in |[A]_q/q|_2, or in the ring setting, |[a]_q / q|^can_2. We take a worst case bound on this value
def worst_case_ATA(n):
    return (2 * n / pi) ** 2 + 1

# this is the amount of required noise to use a HintLWE reduction when the original secret width is sigma, the secret width after hints in sigma_prime.
# t: the number of allowed decryptions
# ATA_bound: a bound on the hint matrix |[A]_q / q|_2^2  
def sigma_i(t, ATA_bound, sigma_prime_sq, sigma_sq):
    sigma_i_sq = log(t, 2) + log(ATA_bound, 2) - log(1/(2 * sigma_prime_sq) - 1/sigma_sq)
    return sigma_i_sq / 2

def hint_lwe_flooding_noise(t, rescaled_noise_bound, target_security, n, sigma_prime_sq, sigma_sq):
    ATA_bound = worst_case_ATA(n)
    sigma_i_ = sigma_i(t, ATA_bound, sigma_prime_sq, sigma_sq)
    rescaled_noise_flooding = bit_security_flooding_noise(rescaled_noise_bound, target_security)
    return max(rescaled_noise_flooding, sigma_i_)

if __name__ == "__main__":
    (logn, logQ, sigma_s, sigma_e) = (13, 216, 3.19, 3.19)
    # this is the security level we find for LWE with these parameters: we take this as the IND-CPA security of the scheme
    original_security = 128.193639623478
    sigma_sq = sigma_s ** 2
    # this is the security level we find for HintLWE with these parameters according to our reduction.
    # if we find a different bit security, we can overwrite this value
    kappa_prime = 112
    # this is the variance of the secret after hints that achieves kappa_prime bits of security: we use this value to calculate the flooding variance
    sigma_prime_sq = 5.05659547963317
    # the tolerable loss vs. the HintLWE security, l in Theorem 10. We use the most efficient reduction, which sets l = 4
    loss = 4
    target_security = kappa_prime - 4 - log(12, 2).n()
    # number of queries to the decryption oracle
    t = 2 ** 6
    max_rescaled_noise_magnitude = 40
    
    n = 2 ** logn
    print(f"# noise flooding to achieve a security level of {target_security}")
    
    # since our results are in terms of |error|_2, we convert this to addition precision loss via sqrt(n) * sigma / |error|_2
    additional_precision_loss = lambda logn, log_sigma, log_total_noise: (logn / 2 + log_sigma - log_total_noise).n()
    
    absolute_noise = {"prior": [], "bit_security": [], "hint_lwe": []}
    
    additional_noise = {"prior": [], "bit_security": [], "hint_lwe": []}
    
    # a (2 norm) bound on the noise due to rescaling, ([A]_q * s - [b]_q ) / q
    rescaling_noise = sqrt(worst_case_ATA(n)) * sqrt(n) * sigma_s + 0.5
    print(f"# cross over point (rescaling noise = rescaled noise) at {log(rescaling_noise, 2).n()}")
    for rescaled_noise_magnitude in range(1, max_rescaled_noise_magnitude + 1):
        # first we calculate a bound on the total noise from the bound on the rescaled noise magnitude.
        total_noise = 2 ** rescaled_noise_magnitude + rescaling_noise
        log_total_noise = log(total_noise, 2)
        
        prior_flooding_noise_ = prior_flooding_noise(t, log_total_noise, target_security)
        absolute_noise["prior"].append((rescaled_noise_magnitude, prior_flooding_noise_))
        prior_precision_loss = additional_precision_loss(logn, prior_flooding_noise_, log_total_noise)
        additional_noise["prior"].append((rescaled_noise_magnitude, prior_precision_loss))
        
        bit_security_flooding_noise_ = bit_security_flooding_noise(log_total_noise, target_security)
        absolute_noise["bit_security"].append((rescaled_noise_magnitude, bit_security_flooding_noise_))
        bit_security_precision_loss = additional_precision_loss(logn, bit_security_flooding_noise_, log_total_noise)
        additional_noise["bit_security"].append((rescaled_noise_magnitude, bit_security_precision_loss))

        hint_lwe_flooding_noise_ = hint_lwe_flooding_noise(t, rescaled_noise_magnitude, target_security, n, sigma_prime_sq, sigma_sq)
        absolute_noise["hint_lwe"].append((rescaled_noise_magnitude, hint_lwe_flooding_noise_))
        hint_lwe_precision_loss = additional_precision_loss(logn, hint_lwe_flooding_noise_, log_total_noise)
        additional_noise["hint_lwe"].append((rescaled_noise_magnitude, hint_lwe_precision_loss))
    
    print("# absolute noise:")
    for key in absolute_noise.keys():
        print(f"#\t{key}:", absolute_noise[key])
        print()
    print()
    print()
    for key in additional_noise.keys():
        print(f"#\t{key}:", additional_noise[key])
        print()
    
## output:
# noise flooding to achieve a security level of 104.415037499279
# cross over point (rescaling noise = rescaled noise) at 20.5220608000984
#       prior: [(1, 78.7295814659711), (2, 78.7295833822018), (3, 78.7295872146556), (4, 78.7295948795326), (5, 78.7296102091645), (6, 78.7296408679396), (7, 78.7297021835352), (8, 78.7298248069093), (9, 78.7300700223942), (10, 78.7305603283614), (11, 78.7315404406811), (12, 78.7334986700261), (13, 78.7374071727539), (14, 78.7451925545611), (15, 78.7606384011157), (16, 78.7910426737623), (17, 78.8499939334301), (18, 78.9611275879052), (19, 79.1606014487134), (20, 79.4920358603678), (21, 79.9882508970452), (22, 80.6500655400265), (23, 81.4456949368854), (24, 82.3315163693188), (25, 83.2708493273730), (26, 84.2395315166613), (27, 85.2236139249774), (28, 86.2155887825687), (29, 87.2115594087853), (30, 88.2095404938323), (31, 89.2085299758863), (32, 90.2080244513626), (33, 91.2077716226587), (34, 92.2076451916894), (35, 93.2075819720496), (36, 94.2075503611909), (37, 95.2075345555017), (38, 96.2075266525922), (39, 97.2075227011212), (40, 98.2075207253817)]

#       bit_security: [(1, 73.7295814659711), (2, 73.7295833822018), (3, 73.7295872146556), (4, 73.7295948795326), (5, 73.7296102091645), (6, 73.7296408679396), (7, 73.7297021835352), (8, 73.7298248069093), (9, 73.7300700223942), (10, 73.7305603283614), (11, 73.7315404406811), (12, 73.7334986700261), (13, 73.7374071727539), (14, 73.7451925545611), (15, 73.7606384011157), (16, 73.7910426737623), (17, 73.8499939334301), (18, 73.9611275879052), (19, 74.1606014487134), (20, 74.4920358603678), (21, 74.9882508970452), (22, 75.6500655400265), (23, 76.4456949368854), (24, 77.3315163693188), (25, 78.2708493273730), (26, 79.2395315166613), (27, 80.2236139249774), (28, 81.2155887825687), (29, 82.2115594087853), (30, 83.2095404938323), (31, 84.2085299758863), (32, 85.2080244513626), (33, 86.2077716226587), (34, 87.2076451916894), (35, 88.2075819720496), (36, 89.2075503611909), (37, 90.2075345555017), (38, 91.2075266525922), (39, 92.2075227011212), (40, 93.2075207253817)]

#       hint_lwe: [(1, 54.2075187496394), (2, 55.2075187496394), (3, 56.2075187496394), (4, 57.2075187496394), (5, 58.2075187496394), (6, 59.2075187496394), (7, 60.2075187496394), (8, 61.2075187496394), (9, 62.2075187496394), (10, 63.2075187496394), (11, 64.2075187496394), (12, 65.2075187496394), (13, 66.2075187496394), (14, 67.2075187496394), (15, 68.2075187496394), (16, 69.2075187496394), (17, 70.2075187496394), (18, 71.2075187496394), (19, 72.2075187496394), (20, 73.2075187496394), (21, 74.2075187496394), (22, 75.2075187496394), (23, 76.2075187496394), (24, 77.2075187496394), (25, 78.2075187496394), (26, 79.2075187496394), (27, 80.2075187496394), (28, 81.2075187496394), (29, 82.2075187496394), (30, 83.2075187496394), (31, 84.2075187496394), (32, 85.2075187496394), (33, 86.2075187496394), (34, 87.2075187496394), (35, 88.2075187496394), (36, 89.2075187496394), (37, 90.2075187496394), (38, 91.2075187496394), (39, 92.2075187496394), (40, 93.2075187496394)]



#       prior: [(1, 64.7075187496394), (2, 64.7075187496394), (3, 64.7075187496394), (4, 64.7075187496394), (5, 64.7075187496394), (6, 64.7075187496394), (7, 64.7075187496394), (8, 64.7075187496394), (9, 64.7075187496394), (10, 64.7075187496394), (11, 64.7075187496394), (12, 64.7075187496394), (13, 64.7075187496394), (14, 64.7075187496394), (15, 64.7075187496394), (16, 64.7075187496394), (17, 64.7075187496394), (18, 64.7075187496394), (19, 64.7075187496394), (20, 64.7075187496394), (21, 64.7075187496394), (22, 64.7075187496394), (23, 64.7075187496394), (24, 64.7075187496394), (25, 64.7075187496394), (26, 64.7075187496394), (27, 64.7075187496394), (28, 64.7075187496394), (29, 64.7075187496394), (30, 64.7075187496394), (31, 64.7075187496394), (32, 64.7075187496394), (33, 64.7075187496394), (34, 64.7075187496394), (35, 64.7075187496394), (36, 64.7075187496394), (37, 64.7075187496394), (38, 64.7075187496394), (39, 64.7075187496394), (40, 64.7075187496394)]

#       bit_security: [(1, 59.7075187496394), (2, 59.7075187496394), (3, 59.7075187496394), (4, 59.7075187496394), (5, 59.7075187496394), (6, 59.7075187496394), (7, 59.7075187496394), (8, 59.7075187496394), (9, 59.7075187496394), (10, 59.7075187496394), (11, 59.7075187496394), (12, 59.7075187496394), (13, 59.7075187496394), (14, 59.7075187496394), (15, 59.7075187496394), (16, 59.7075187496394), (17, 59.7075187496394), (18, 59.7075187496394), (19, 59.7075187496394), (20, 59.7075187496394), (21, 59.7075187496394), (22, 59.7075187496394), (23, 59.7075187496394), (24, 59.7075187496394), (25, 59.7075187496394), (26, 59.7075187496394), (27, 59.7075187496394), (28, 59.7075187496394), (29, 59.7075187496394), (30, 59.7075187496394), (31, 59.7075187496394), (32, 59.7075187496394), (33, 59.7075187496394), (34, 59.7075187496394), (35, 59.7075187496394), (36, 59.7075187496394), (37, 59.7075187496394), (38, 59.7075187496394), (39, 59.7075187496394), (40, 59.7075187496394)]

#       hint_lwe: [(1, 40.1854560333077), (2, 41.1854541170770), (3, 42.1854502846233), (4, 43.1854426197462), (5, 44.1854272901144), (6, 45.1853966313393), (7, 46.1853353157436), (8, 47.1852126923696), (9, 48.1849674768846), (10, 49.1844771709174), (11, 50.1834970585977), (12, 51.1815388292528), (13, 52.1776303265250), (14, 53.1698449447178), (15, 54.1543990981632), (16, 55.1239948255165), (17, 56.0650435658488), (18, 56.9539099113736), (19, 57.7544360505654), (20, 58.4230016389110), (21, 58.9267866022336), (22, 59.2649719592524), (23, 59.4693425623934), (24, 59.5835211299600), (25, 59.6441881719059), (26, 59.6755059826175), (27, 59.6914235743014), (28, 59.6994487167101), (29, 59.7034780904935), (30, 59.7054970054465), (31, 59.7065075233925), (32, 59.7070130479163), (33, 59.7072658766202), (34, 59.7073923075894), (35, 59.7074555272292), (36, 59.7074871380880), (37, 59.7075029437771), (38, 59.7075108466866), (39, 59.7075147981576), (40, 59.7075167738972)]

        
        
        
        
    
    

    
    




    

