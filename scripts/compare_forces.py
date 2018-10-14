#!/anaconda/bin/python
# Compare two files with forces

import sys, math
import numpy as np

# Input: files with time series of forces
argv = sys.argv
argc = len(argv)
if ( argc < 3 ):
  print("Usage: python compare_forces.py <file1> <file2>")
  exit(1)
file_1 = argv[1]
file_2 = argv[2]

# Average magnitude of force components
avg_abs_f1 = np.zeros(3, dtype=float)
avg_abs_f2 = np.zeros(3, dtype=float)
max_abs_f1 = np.zeros(3, dtype=float)
max_abs_f2 = np.zeros(3, dtype=float)
num_samples_f1 = np.zeros(3, dtype=int)
num_samples_f2 = np.zeros(3, dtype=int)

# Largest error in a force/virial component
max_abs_df  = 0.0
max_abs_dxi = 0.0

avg_df        = 0.0
avg_df_sq     = 0.0
avg_abs_df    = 0.0
num_samples_df = 0

line_number = 0

with open(file_1) as fp_1, open(file_2) as fp_2:
	# Read lines
	for line_1, line_2 in zip(fp_1, fp_2):

#	line_1 = fp_1.readline()
#	line_2 = fp_2.readline()
		line_number += 1

		# Tokenize lines
		tokens_1 = line_1.strip().split()
		tokens_2 = line_2.strip().split()
		num_tokens = len(tokens_1)
		if ( num_tokens != len(tokens_2) ):
			print("ERROR: Different numbers of tokens at line {}:\n{}\n{}\n".format( \
						line_number, line_1, line_2))
			exit(1)

		# Types of lines
		if ( num_tokens == 0 or tokens_1[0] == '#' or \
				 (tokens_1[0] == 't=' and tokens_2[0]  == 't=') or \
				 (tokens_1[0]  == 'box=' and tokens_2[0]  == 'box=') or \
				 (tokens_1[0]  == 'N_total=' and tokens_2[0]  == 'N_total=') ):
			continue
		elif ( tokens_1[0] == 'virial=' and tokens_2[0] == 'virial=' ):
			# Virial tensor
			for i in range(1, num_tokens):
				abs_dxi = abs( float(tokens_1[i]) - float(tokens_2[i]) )
				if ( abs_dxi > max_abs_dxi ):
					max_abs_dxi = abs_dxi
		elif ( int(tokens_1[0]) == int(tokens_2[0]) ):
			# Force
			for i in range(1, num_tokens):
				f_1 = float(tokens_1[i])
				f_2 = float(tokens_2[i])
				# Difference
				df = f_1 - f_2
				abs_df = abs(df)
				# Maximum
				if ( abs_df > max_abs_df ):
					max_abs_df = abs_df
				# Averages
				if ( f_1 != 0.0 or f_2 != 0.0 ):
					avg_df        += df
					avg_df_sq     += df*df
					avg_abs_df    += abs_df
					num_samples_df += 1

				# Statistics of nonzero components
				if ( f_1 != 0.0 ):
					abs_f1 = abs(f_1)

					if ( abs_f1 > max_abs_f1[i-1] ):
						max_abs_f1[i-1] = abs_f1

					avg_abs_f1[i-1] += abs_f1
					num_samples_f1[i-1] += 1
				if ( f_2 != 0.0 ):
					abs_f2 = abs(f_2)

					if ( abs_f2 > max_abs_f2[i-1] ):
						max_abs_f2[i-1] = abs_f2

					avg_abs_f2[i-1] += abs_f2
					num_samples_f2[i-1] += 1
		else:
			print("ERROR: Could not understand the content at line {}:\n{}\n{}\n".format( \
						line_number, line_1, line_2))
			exit(1)

# Finish statistics
avg_abs_f1 /= num_samples_f1.astype(float)
avg_abs_f2 /= num_samples_f2.astype(float)

avg_df     /= float(num_samples_df)
avg_df_sq  /= float(num_samples_df)
avg_abs_df /= float(num_samples_df)

var_df = float(num_samples_df)/(float(num_samples_df)+1)*(avg_df_sq - avg_df*avg_df)
std_dev_df = math.sqrt(var_df)

var_abs_df = float(num_samples_df)/(float(num_samples_df)+1)*(avg_df_sq - avg_abs_df*avg_abs_df)
std_dev_abs_df = math.sqrt(var_abs_df)

# Output
print("Compared forces")

dim_names = np.array(['x', 'y', 'z'], dtype=str)
dim = dim_names.shape[0]
print("    avg(|f_1,d|)  avg(|f_2,d|)  ")
for d in range(0,dim):
	print("{}:  {}  {}".format( dim_names[d], avg_abs_f1[d], avg_abs_f2[d] ))

print("    max(|f_1,d|)  max(|f_2,d|)  ")
for d in range(0,dim):
	print("{}:  {}  {}".format( dim_names[d], max_abs_f1[d], max_abs_f2[d] ))

print("delta_f = f_i(1) - f_i(2)   [for f(1) and f(2) not both zero]")
print("  delta_f")
print("    avg     = {:.3e}".format(avg_df))
print("    std_dev = {:.3e}".format(std_dev_df))
print("  abs(delta_f)")
print("    avg     = {:.3e}".format(avg_abs_df))
print("    std_dev = {:.3e}".format(std_dev_abs_df))
print("    max     = {:.3e}".format(max_abs_df))
print("")
print("delta_Xi = Xi_ij(1) - Xi_ij(2)   [for Xi_ij(1) and Xi_ij(2) not both zero]")
print("  abs(delta_Xi)")
print("    max = {:.3e}".format(max_abs_dxi))
