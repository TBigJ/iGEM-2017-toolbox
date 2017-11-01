import numpy as np
from numpy import log
### Parameters for the simulation

# All time units are in seconds,
# all concentration units micromols

# Literature parameters; see parameter_set.xlsx for source
H_lit     = 2.0/60
K_lit     = 3.7 
P_tot_lit = 0.1 # which is with in range: 0.075 to 0.25
k_lit = 11.3 / 60 # this is called V_{max} in the literature

# delta is computed based on ecoli doubling time of 20 minutes from
# http://sandwalk.blogspot.com/2008/05/dna-replication-in-e-coli-solution.html
# we know that ln(2)/delta = [doubling time]
# and hence delta = ln(2)/doubling time, which is 1200 seconds = 20 mins
delta_lit = log(2)/1200

def literature_parameter_set(n_proteins):
  return {
    'H': np.repeat(H_lit, n_proteins),
    'K': np.repeat(K_lit, n_proteins),
    'P_tot': P_tot_lit,
    'k': k_lit,
    'delta': delta_lit,
  }