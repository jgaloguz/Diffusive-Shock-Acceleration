# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import sys

# Load likelihood and diffusion and compute artificial drift
h0 = np.loadtxt("dsa_results/likelihood_tidx_0.dat")
X0 = h0[:,0]
L0 = h0[:,1]
K0 = h0[:,2]
dX0 = X0[1] - X0[0]
dL0 = (np.log(L0)[2:] - np.log(L0)[:-2]) / (2.0 * dX0)
h1 = np.loadtxt("dsa_results/likelihood_tidx_1.dat")
X1 = h1[:,0]
L1 = h1[:,1]
K1 = h1[:,2]
dX1 = X1[1] - X1[0]
dL1 = (np.log(L1)[2:] - np.log(L1)[:-2]) / (2.0 * dX1)
h2 = np.loadtxt("dsa_results/likelihood_tidx_2.dat")
X2 = h2[:,0]
L2 = h2[:,1]
K2 = h2[:,2]
dX2 = X2[1] - X2[0]
dL2 = (np.log(L2)[2:] - np.log(L2)[:-2]) / (2.0 * dX2)
h3 = np.loadtxt("dsa_results/likelihood_tidx_3.dat")
X3 = h3[:,0]
L3 = h3[:,1]
K3 = h3[:,2]
dX3 = X3[1] - X3[0]
dL3 = (np.log(L3)[2:] - np.log(L3)[:-2]) / (2.0 * dX3)

# Plot data
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')
ax.semilogy(X0, L0, linewidth=3, label="$\\tau$ = $t_0$")
ax.semilogy(X1, L1, linewidth=3, label="$\\tau$ = $t_1$")
ax.semilogy(X2, L2, linewidth=3, label="$\\tau$ = $t_2$")
ax.set_xlabel('$x$ (au)', fontsize=20)
ax.set_ylabel('$h(x,0)$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

plt.savefig("likelihood_sol_slices.png")
plt.show()
plt.close(fig)

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')
ax.plot(X0[1:-1], 2.0 * K0[1:-1] * dL0, linewidth=3, label="$\\tau$ = $t_0$")
ax.plot(X1[1:-1], 2.0 * K1[1:-1] * dL1, linewidth=3, label="$\\tau$ = $t_1$")
ax.plot(X2[1:-1], 2.0 * K2[1:-1] * dL2, linewidth=3, label="$\\tau$ = $t_2$")
ax.set_xlabel('$x$ (au)', fontsize=20)
ax.set_ylabel('$2 \\kappa(x) \\partial_x \\log h(x, 0)$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

plt.savefig("likelihood_drift_slices.png")
plt.show()
plt.close(fig)