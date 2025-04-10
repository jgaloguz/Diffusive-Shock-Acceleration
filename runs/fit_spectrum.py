# Script to fit momentum spectrum with an exponential curve (linear in log-log view)

# Import libraries
import numpy as np
import sys

# Import data
data = np.loadtxt(sys.argv[1])

# Log momentum and flux
p = data[:,0]
J = data[:,1]

# Remove zero values
zero_indices = np.where(J == 0.0)
logp = np.log(np.delete(p, zero_indices))
logJ = np.log(np.delete(J, zero_indices))

# Fit with a line
par, cov = np.polyfit(logp, logJ, 1, cov=True)

# Output slope and variance of slope
print("slope =", par[0])
print("std =", np.sqrt(cov[0][0]))
