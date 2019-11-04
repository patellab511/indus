#!/usr/bin/env python
# Compare two files with forces (Indus code format)
# - Note: numerical derivatives and analytic derivatives produce files
#         with the same format

import sys, math
import numpy as np

# Plotting with matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('text', **{'usetex': True})
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})  # LaTeX font
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # AMS Math (includes \text{})

###

# Input: files with time series of forces
argv = sys.argv
argc = len(argv)
if ( argc < 3 ):
  print("Usage: python compare_forces.py <file1> <file2>")
  exit(1)
file_1 = argv[1]
file_2 = argv[2]

dim = 3

# Samples for virial
xi_1_samples = []
xi_2_samples = []
xi_1 = np.zeros((dim,dim), dtype=np.float)
xi_2 = np.zeros((dim,dim), dtype=np.float)

# Samples for forces
f_1_samples = []
f_2_samples = []
f_1    = np.zeros(dim, dtype=np.float)
f_2    = np.zeros(dim, dtype=np.float)
f_zero = np.zeros(dim, dtype=np.float)

line_number = 0

tokens_to_skip = [ "t=", "box=", "N_total=" ]

# Token offsets for reading virial/forces
xi_offset = 1
f_offset = 1

with open(file_1) as fp_1, open(file_2) as fp_2:
	# Read lines
	for line_1, line_2 in zip(fp_1, fp_2):
		line_number += 1

		# Tokenize lines
		tokens_1 = line_1.strip().split()
		tokens_2 = line_2.strip().split()
		num_tokens = len(tokens_1)
		if ( num_tokens != len(tokens_2) ):
			raise RuntimeError("ERROR: Different numbers of tokens at line {}:\n{}\n{}\n".format( \
			                   line_number, line_1, line_2))

		# Types of lines
		if ( num_tokens == 0 or tokens_1[0] == '#'
				 or ( tokens_1[0] in tokens_to_skip ) 
				 or ( tokens_2[0] in tokens_to_skip )
		):
			# Skip this line
			continue
		elif ( tokens_1[0] == 'virial=' and tokens_2[0] == 'virial=' ):
			# Virial tensor
			for a in range(0,dim):
				for b in range(0,dim):
					i = a*dim + b + xi_offset
					xi_1[a][b] = float(tokens_1[i])
					xi_2[a][b] = float(tokens_2[i])
			xi_1_samples.append( xi_1.copy() )
			xi_2_samples.append( xi_2.copy() )
		elif ( int(tokens_1[0]) == int(tokens_2[0]) ):
			for d in range(0,dim):
				i = d + f_offset
				f_1[d] = float(tokens_1[i])
				f_2[d] = float(tokens_2[i])

			if ( (not np.array_equal(f_1, f_zero)) or (not np.array_equal(f_2, f_zero)) ):
				# Save samples
				f_1_samples.append( f_1.copy() )
				f_2_samples.append( f_2.copy() )
		else:
			print("ERROR: Could not understand the content at line {}:\n{}\n{}\n".format( \
						line_number, line_1, line_2))
			exit(1)

# Convert to numpy arrays
f_1_samples = np.array(f_1_samples)
f_2_samples = np.array(f_2_samples)
xi_1_samples = np.array(xi_1_samples)
xi_2_samples = np.array(xi_2_samples)

# Element-wise absolute value of forces
abs_f_1 = np.abs( f_1_samples )
abs_f_2 = np.abs( f_2_samples )
avg_abs_f_1 = np.average(abs_f_1, axis=0)
avg_abs_f_2 = np.average(abs_f_2, axis=0)
max_abs_f_1 = np.max(abs_f_1, axis=0)
max_abs_f_2 = np.max(abs_f_2, axis=0)

# Difference between forces (and its absolute value)
df     = f_1_samples - f_2_samples
abs_df = np.abs(df)
avg_df     = np.average(df)
std_dev_df = np.std(df, ddof=1)
avg_abs_df     = np.average(abs_df)
std_dev_abs_df = np.std(abs_df, ddof=1)

indices_max_df = np.unravel_index(np.argmax(abs_df, axis=None), abs_df.shape)
max_abs_df     = abs_df[indices_max_df]

# Statistics of delta_xi
delta_xi = xi_1_samples - xi_2_samples
abs_delta_xi = np.abs(delta_xi)
max_abs_delta_xi = np.max(abs_delta_xi) 

# Output
print("Compared forces")

axis_names = np.array(['x', 'y', 'z'], dtype=str)
print("    avg(|f_1,d|)  avg(|f_2,d|)  ")
for d in range(0,dim):
	print("{}:  {}  {}".format( axis_names[d], avg_abs_f_1[d], avg_abs_f_2[d] ))

print("    max(|f_1,d|)  max(|f_2,d|)  ")
for d in range(0,dim):
	print("{}:  {}  {}".format( axis_names[d], max_abs_f_1[d], max_abs_f_2[d] ))

print("delta_f = f_i(1) - f_i(2)   [for f(1) and f(2) not both zero]")
print("  delta_f")
print("    avg     = {:.3e}".format(avg_df))
print("    std_dev = {:.3e}".format(std_dev_df))
print("  abs(delta_f)")
print("    avg     = {:.3e}".format(avg_abs_df))
print("    std_dev = {:.3e}".format(std_dev_abs_df))
print("    max     = {:.3e} (f_1 = {:.3e}, f_2 = {:.3e})".format(
                     max_abs_df, f_1_samples[indices_max_df], f_2_samples[indices_max_df]))
print("")
print("delta_Xi = Xi_ij(1) - Xi_ij(2)   [for Xi_ij(1) and Xi_ij(2) not both zero]")
print("  abs(delta_Xi)")
print("    max = {:.3e}".format(max_abs_delta_xi))

make_plots = True
if ( make_plots ):
	# Parity plot: f1 vs. f2 for each axis
	for d in range(0, dim):
		plt.clf()
		fig, ax = plt.subplots()

		ax.plot( f_1_samples[:,d], f_2_samples[:,d], '.', markersize=3, label="Data" )
		x1,x2 = ax.get_xlim()
		y1,y2 = ax.get_ylim()
		min_xy = min( x1, x2, y1, y2 )
		max_xy = max( x1, x2, y1, y2 )

		parity = np.linspace(min_xy, max_xy)
		ax.plot( parity, parity, '-', color='black', label='$y = x$', 
		         linestyle='dashed', linewidth=1 )

		ax.set_xlabel("$F_{" + axis_names[d] + ",1}$ (kJ/mol$\\cdot$nm)")
		ax.set_ylabel("$F_{" + axis_names[d] + ",2}$ (kJ/mol$\\cdot$nm)")
		plt.axis('scaled')
		plt.legend()
		plt.tight_layout()
		fig.savefig("plot_parity_force_" + axis_names[d] + ".png", dpi=300)
		plt.close(fig)
