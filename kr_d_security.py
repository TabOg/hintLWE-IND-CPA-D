# this code corresponds to Section 6.3
from sage.all import log, e, sqrt
from ind_cpa_d_security import worst_case_ATA, sigma_i

def bit_security_flooding_noise(t, logE, original_security, target_security):
    l = original_security - target_security
    assert(l >= 1)
    loss = original_security - target_security
    total_flooding_noise = log(2 * original_security * t * log(e, 2), 2) / 2 + logE - log(l, 2)
    return total_flooding_noise.n()


def hint_lwe_flooding_noise(t, rescaled_noise_bound, target_security, n, sigma_prime_sq, sigma_sq, loss=1):
    ATA_bound = worst_case_ATA(n)
    sigma_i_ = sigma_i(t, ATA_bound, sigma_prime_sq, sigma_sq)
    
    rescaled_noise_flooding = bit_security_flooding_noise(t, rescaled_noise_bound, target_security + loss, target_security)
    
    return max(rescaled_noise_flooding, sigma_i_)

if __name__ == "__main__":
    (logn, logQ, sigma_s, sigma_e) = (13, 216, 3.19, 3.19)
    # this is the security level we find for LWE with these parameters: we take this as the IND-CPA security of the scheme
    original_security = 128.193639623478
    sigma_sq = sigma_s ** 2
    # this is the security level we find for HintLWE with these parameters according to our reduction.
    # if we find a different bit security, we can overwrite this value
    kappa_prime = 126.198049311874
    # this is the variance of the secret after hints that achieves kappa_prime bits of security: we use this value to calculate the flooding variance
    sigma_prime_sq = 1.17787445304450
    # the tolerable loss vs. the HintLWE security, l in Theorem 11. We use the most efficient rediction, which sets l = 1
    loss = 1
    target_security = kappa_prime - loss
    # number of queries to the decryption oracle
    t = 2 ** 6
    max_rescaled_noise_magnitude = 40
    
    n = 2 ** logn
    print(f"# noise flooding to achieve a security level of {target_security}")
    # since our results are in terms of |error|_2, we convert this to addition precision loss via sqrt(n) * sigma / |error|_2
    additional_precision_loss = lambda logn, log_sigma, log_total_noise: (logn / 2 + log_sigma - log_total_noise).n()
    
    absolute_noise = {"bit_security": [], "hint_lwe": []}
    additional_noise = {"bit_security": [], "hint_lwe": []}
    
    # a (2 norm) bound on the noise due to rescaling, ([A]_q * s - [b]_q ) / q
    rescaling_noise = sqrt(worst_case_ATA(n)) * sqrt(n) * sigma_s + 0.5
    print(f"# cross over point (rescaling noise = rescaled noise) at {log(rescaling_noise, 2).n()}")
    
    for rescaled_noise_magnitude in range(1, max_rescaled_noise_magnitude + 1):
        # first we calculate a bound on the total noise from the bound on the rescaled noise magnitude. All in 2 norm
        total_noise = 2 ** rescaled_noise_magnitude + rescaling_noise
        log_total_noise = log(total_noise, 2)
        
        bit_security_flooding_noise_ = bit_security_flooding_noise(t, log_total_noise, original_security, target_security)
        absolute_noise["bit_security"].append((rescaled_noise_magnitude, bit_security_flooding_noise_))
        bit_security_precision_loss = additional_precision_loss(logn, bit_security_flooding_noise_, log_total_noise)
        additional_noise["bit_security"].append((rescaled_noise_magnitude, bit_security_precision_loss))
        
        hint_lwe_flooding_noise_ = hint_lwe_flooding_noise(t, rescaled_noise_magnitude, target_security, n, sigma_prime_sq, sigma_sq, loss=loss)
        absolute_noise["hint_lwe"].append((rescaled_noise_magnitude, hint_lwe_flooding_noise_.n()))
        hint_lwe_precision_loss = additional_precision_loss(logn, hint_lwe_flooding_noise_, log_total_noise)
        additional_noise["hint_lwe"].append((rescaled_noise_magnitude, hint_lwe_precision_loss.n()))
        
    print("# absolute noise:")
    for key in absolute_noise.keys():
        print(f"#\t{key}:", absolute_noise[key])
        print()
    print()
    print()
    for key in additional_noise.keys():
        print(f"#\t{key}:", additional_noise[key])
        print()
# output:
# noise flooding to achieve a security level of 125.198049311874
# cross over point (rescaling noise = rescaled noise) at 20.5220608000984
# absolute noise:
#       bit_security: [(1, 26.2046960108318), (2, 26.2046979270625), (3, 26.2047017595163), (4, 26.2047094243933), (5, 26.2047247540252), (6, 26.2047554128003), (7, 26.2048167283960), (8, 26.2049393517700), (9, 26.2051845672549), (10, 26.2056748732221), (11, 26.2066549855419), (12, 26.2086132148868), (13, 26.2125217176146), (14, 26.2203070994218), (15, 26.2357529459764), (16, 26.2661572186231), (17, 26.3251084782908), (18, 26.4362421327659), (19, 26.6357159935742), (20, 26.9671504052285), (21, 27.4633654419060), (22, 28.1251800848872), (23, 28.9208094817461), (24, 29.8066309141795), (25, 30.7459638722337), (26, 31.7146460615220), (27, 32.6987284698381), (28, 33.6907033274294), (29, 34.6866739536461), (30, 35.6846550386930), (31, 36.6836445207471), (32, 37.6831389962233), (33, 38.6828861675194), (34, 39.6827597365502), (35, 40.6826965169104), (36, 41.6826649060516), (37, 42.6826491003624), (38, 43.6826411974529), (39, 44.6826372459820), (40, 45.6826352702424)]

#       hint_lwe: [(1, 15.9085894304852), (2, 15.9085894304852), (3, 15.9085894304852), (4, 15.9085894304852), (5, 15.9085894304852), (6, 15.9085894304852), (7, 15.9085894304852), (8, 15.9085894304852), (9, 16.2541560864811), (10, 17.2541560864811), (11, 18.2541560864811), (12, 19.2541560864811), (13, 20.2541560864811), (14, 21.2541560864811), (15, 22.2541560864811), (16, 23.2541560864811), (17, 24.2541560864811), (18, 25.2541560864811), (19, 26.2541560864811), (20, 27.2541560864811), (21, 28.2541560864811), (22, 29.2541560864811), (23, 30.2541560864811), (24, 31.2541560864811), (25, 32.2541560864811), (26, 33.2541560864811), (27, 34.2541560864811), (28, 35.2541560864811), (29, 36.2541560864811), (30, 37.2541560864811), (31, 38.2541560864811), (32, 39.2541560864811), (33, 40.2541560864811), (34, 41.2541560864811), (35, 42.2541560864811), (36, 43.2541560864811), (37, 44.2541560864811), (38, 45.2541560864811), (39, 46.2541560864811), (40, 47.2541560864811)]



#       bit_security: [(1, 12.1826332945001), (2, 12.1826332945001), (3, 12.1826332945001), (4, 12.1826332945001), (5, 12.1826332945001), (6, 12.1826332945001), (7, 12.1826332945001), (8, 12.1826332945001), (9, 12.1826332945001), (10, 12.1826332945001), (11, 12.1826332945001), (12, 12.1826332945001), (13, 12.1826332945001), (14, 12.1826332945001), (15, 12.1826332945001), (16, 12.1826332945001), (17, 12.1826332945001), (18, 12.1826332945001), (19, 12.1826332945001), (20, 12.1826332945001), (21, 12.1826332945001), (22, 12.1826332945001), (23, 12.1826332945001), (24, 12.1826332945001), (25, 12.1826332945001), (26, 12.1826332945001), (27, 12.1826332945001), (28, 12.1826332945001), (29, 12.1826332945001), (30, 12.1826332945001), (31, 12.1826332945001), (32, 12.1826332945001), (33, 12.1826332945001), (34, 12.1826332945001), (35, 12.1826332945001), (36, 12.1826332945001), (37, 12.1826332945001), (38, 12.1826332945001), (39, 12.1826332945001), (40, 12.1826332945001)]

#       hint_lwe: [(1, 1.88652671415356), (2, 1.88652479792285), (3, 1.88652096546906), (4, 1.88651330059205), (5, 1.88649797096017), (6, 1.88646731218507), (7, 1.88640599658941), (8, 1.88628337321538), (9, 2.23160481372627), (10, 3.23111450775908), (11, 4.23013439543935), (12, 5.22817616609443), (13, 6.22426766336661), (14, 7.21648228155941), (15, 8.20103643500482), (16, 9.17063216235817), (17, 10.1116809026904), (18, 11.0005472482153), (19, 11.8010733874071), (20, 12.4696389757527), (21, 12.9734239390752), (22, 13.3116092960940), (23, 13.5159798992351), (24, 13.6301584668017), (25, 13.6908255087475), (26, 13.7221433194592), (27, 13.7380609111431), (28, 13.7460860535518), (29, 13.7501154273351), (30, 13.7521343422882), (31, 13.7531448602341), (32, 13.7536503847579), (33, 13.7539032134618), (34, 13.7540296444311), (35, 13.7540928640708), (36, 13.7541244749296), (37, 13.7541402806188), (38, 13.7541481835283), (39, 13.7541521349993), (40, 13.7541541107388)]
    