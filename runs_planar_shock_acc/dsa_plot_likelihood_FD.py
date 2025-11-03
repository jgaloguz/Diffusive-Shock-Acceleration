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
dX = X[1:] - X[:-1]
Cl = 1.0 / (dX[1:] + dX[:-1]) - 1.0 / dX[:-1]
Cc = 1.0 / dX[:-1] - 1.0 / dX[1:]
Cr = 1.0 / dX[1:] - 1.0 / (dX[1:] + dX[:-1])
U = coeffs[:,1]
K = coeffs[:,2]
dU = coeffs[:,3]

# Plot solution
t_idxs = [1, 70, 90, 97, 100]
t = [Tf_FD * (Nt_out - 1 - idx) / (Nt_out - 1) for idx in t_idxs]
h = np.loadtxt("dsa_results/likelihood_solution_FD.dat")

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for idx in range(len(t_idxs)):
   ax.semilogy(X, h[t_idxs[idx],:], linewidth=3, label="$t=${:.2f}".format(t[idx]))
ax.set_xlabel('$x$ (au)', fontsize=20)
ax.set_ylabel('$h(x,t)$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlim(-12.0,3.0)
ax.legend(fontsize=20)

plt.savefig("likelihood_sol_slices.png")
plt.show()
plt.close(fig)

# Plot drift
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for idx in range(len(t_idxs)):
   logh = np.log(h[t_idxs[idx],:])
   dlogh = Cl * logh[:-2] + Cc * logh[1:-1] + Cr * logh[2:]
   ax.plot(X[1:-1], 2.0 * K[1:-1] * dlogh, linewidth=3, label="$t=${:.2f}".format(t[idx]))
ax.set_xlabel('$x$ (au)', fontsize=20)
ax.set_ylabel('$2 \\kappa(x) \\partial_x \\log h(x,t)$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlim(-12.0,3.0)
ax.legend(fontsize=20)

plt.savefig("likelihood_drift_slices.png")
plt.show()
plt.close(fig)