# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import sys

# Load parameters
params = np.loadtxt("dsa_results/likelihood_params_FD.dat")
Nx_FD = params[0]
Nx_lft = params[1]
Nx_ctr = params[2]
Nx_rgt = params[3]
X1_FD = params[4]
X2_FD = params[5]
X3_FD = params[6]
X4_FD = params[7]
dt_FD = params[8]
Tf_FD = params[9]
Nt_out = params[10]
alpha = params[11]

# Load coefficients
coeffs = np.loadtxt("dsa_results/likelihood_coeffs_FD.dat")
X = coeffs[:,0]
U = coeffs[:,1]
K = coeffs[:,2]
dU = coeffs[:,3]

# Plot solution
h = np.loadtxt("dsa_results/likelihood_solution_FD.dat")
t_idxs = [0, 3, 10, 30, 100]

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for idx in range(len(t_idxs)):
   ax.semilogy(X, h[t_idxs[idx],:], linewidth=3)
ax.set_xlabel('$x$ (au)', fontsize=20)
ax.set_ylabel('$h(x,t)$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

plt.savefig("likelihood_sol_slices.png")
plt.show()
plt.close(fig)
