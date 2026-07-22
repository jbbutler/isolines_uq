# Bias Correction Simulation Nomenclature

<bias correction strategy>_<tube construction>_<sample splitting strategy>.R

bias correction strategy:
1. oracleregion: compute conservative region to search for true isoline, using the true relative bias
2. oracleextrapregion: compute conservative region to search for true isoline, using the true relative bias but at a point further from the tail where the bias is assumed to be strictly worse
3. estextrapregion: compute conservative region to search for true isoline, using the estimated relative bias but at a point further from the tail where the bias is assumed to be strictly worse

tube construction:
1. asymmetric: estimating c minus alpha and c plus alpha
2. symmetric: estimating a single c alpha

sample splitting strategy (refer to overleaf):