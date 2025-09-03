# IND-CPA-D and KR-D Security with Reduced Noise from the HintLWE Problem
This repository contains code to generate the data from the above paper.
## Dependencies
This repository uses [sage](https://www.sagemath.org/). In addition, we include a specific commit of the [lattice estimator](https://github.com/malb/lattice-estimator/tree/8f1ff7e20a4d3391e3badff1d76825314db225bc) as a submodule: after cloning this repository, this can be downloaded via
```
git submodule update --init --recursive
```
Depending on your environment, it may be necessary to add the path to the estimator to your python path. We achieve this with
```
export PYTHONPATH="$PYTHONPATH:$(pwd)/lattice_estimator"
```
## Files
We include three programs, corresponding to Figures 1, 2, and 3 respectively. These are given by:
- [hintLWE_security.py](hintLWE_security.py) Fig 1. Calculates the bit security of the HintLWE problem according to our reduction from LWE. 
- [ind_cpa_d_security.py](ind_cpa_d_security.py) Fig 2. Calculates the noise necessary to achieve IND-CPA-D security under both a pure bit security and HintLWE approach. This file uses a security estimate and required flooding noise for HintLWE (Decision): we generate these with [hintLWE_security.py](hintLWE_security.py).
- [kr_d_security](kr_d_security.py) Fig 3. Calculates the noise necessary to achieve KR-D security under both a pure bit security and HintLWE approach. Again this file uses a security estimate and required flooding noise for HintLWE (Search): we generate these with [hintLWE_security.py](hintLWE_security.py).

We run these files with python, using e.g.
```
python3 hintLWE_security.py 
```

## See Also
We target parameters and security model used in the working version of the [HE Community Standard](https://eprint.iacr.org/2024/463): further details can be found in the [corresponding repository](https://github.com/gong-cr/FHE-Security-Guidelines/).

