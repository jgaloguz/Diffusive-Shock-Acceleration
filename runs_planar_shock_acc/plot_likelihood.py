# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import data
likelihood = np.loadtxt("dsa_results/likelihood_expectation.dat")
Nx = np.size(likelihood, 0)
Ny = np.size(likelihood, 1)
X = np.arange(Nx)
Y = np.arange(Ny)
XX, YY = np.meshgrid(X, Y)

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

ax.pcolormesh(XX, YY, likelihood)
ax.set_xlabel('x', fontsize=20)
ax.tick_params(axis='x', labelsize=16)
ax.set_ylabel('p', fontsize=20)
ax.tick_params(axis='y', labelsize=16)

plt.savefig("dsa_likelihood_expectation.png")
plt.show()
plt.close(fig)