#!/usr/bin/env python3

import os, sys, shutil, subprocess
import numpy as np

# Use a shell script to read the values
read_virial_script = "read_virial.sh"
shell_command = ["bash", read_virial_script]
proc = subprocess.Popen(shell_command)
proc.wait()

# Load values
indus_virial  = np.loadtxt("indus_virial.out", comments='#', dtype=np.float)
plumed_virial = 0.5 * np.loadtxt("plumed_2x_virial.out", comments='#', dtype=np.float)

# Compute percent error
err     = (plumed_virial/indus_virial - 1.0) * 100.0
abs_err = np.abs(err)

avg_err = err.mean()
std_err = err.std()

avg_abs_err = abs_err.mean()
std_abs_err = abs_err.std()

print("Error (%)")
print("  Relative: {} +/- {}".format(avg_err, std_err))
print("  Absolute: {} +/- {}".format(avg_abs_err, std_abs_err))
