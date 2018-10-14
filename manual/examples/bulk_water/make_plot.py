#!/anaconda/bin/python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from math  import sqrt
from matplotlib import cm

# matplotlib settings for LaTeX
plt.rc('text', **{'usetex': True})
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})  # LaTeX font
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # AMS Math (includes \text{})

#################
### Load data ###
#################

wham_data_file = "F_Ntilde_WHAM.out"
wham_data = np.loadtxt( wham_data_file, comments='#', dtype=float )
ntilde_v = wham_data[:,0]
f_wham   = wham_data[:,1]


########################
### Plot F(Ntilde_v) ###
########################

plt.clf()
fig, ax = plt.subplots()

mask_f_finite = np.isfinite( f_wham )
min_f = np.min( f_wham[mask_f_finite] )

ax.plot( ntilde_v, f_wham, \
         marker='o', markersize='1', linestyle='none', color='black' )

# x/y-axis bounds
#ax.set_xlim( 0.0, max_ntilde_nn )
#ax.set_ylim( 0.0, max_ntilde_E )

# Labels
ax.set_xlabel("$ \\tilde{N} $", fontsize='x-large')
ax.set_ylabel("$ \\beta F_v = - \\ln P_v( \\tilde{N} ) $", fontsize='x-large')

# Save 
fig.set_size_inches(3.5, 3.0)
plt.tight_layout()
fig.savefig( 'plot_F_v_Ntilde_WHAM.eps', format='eps' )
