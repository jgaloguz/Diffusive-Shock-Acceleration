# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import data
Nt = 5
tf = 5.77548
dsa_analytic_pos = np.loadtxt("dsa_results/dsa_analytic_pos_{:d}.dat".format(Nt-1))
dsa_analytic_mom = np.loadtxt("dsa_results/dsa_analytic_mom_{:d}.dat".format(Nt-1))
# which = sys.argv[1]
# dsa_simulated_pos = [np.loadtxt("dsa_results/dsa_" + which + "_pos_{:d}_pp.dat".format(t)) for t in range(Nt)]
# dsa_simulated_mom = [np.loadtxt("dsa_results/dsa_" + which + "_mom_{:d}_pp.dat".format(t)) for t in range(Nt)]

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax1 = fig.add_subplot(211, projection='rectilinear')

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
ax1.semilogy(dsa_analytic_pos[:,0], dsa_analytic_pos[:,1],
             color=colors[Nt-1], label="t = {:.2e}".format(tf))
# for t in range(Nt):
#    ax1.semilogy(dsa_simulated_pos[t][:,0], dsa_simulated_pos[t][:,1],
#                 color=colors[t], linestyle="", marker="o")
ax1.set_xlabel('x', fontsize=20)
ax1.set_ylabel('N', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
# ax1.set_xlim(79.0, 84.0)
# ax1.set_ylim(1.0e-3, 1.5e0)
ax1.legend(fontsize=20)

ax2 = fig.add_subplot(212, projection='rectilinear')

ax2.loglog(dsa_analytic_mom[:,0], dsa_analytic_mom[:,1],
           color=colors[Nt-1], label="t = {:.2e}".format(tf))
# for t in range(Nt):
#    ax2.loglog(dsa_simulated_mom[t][:,0], dsa_simulated_mom[t][:,1],
#               color=colors[t], linestyle="", marker="s")
ax2.set_xlabel('p', fontsize=20)
ax2.set_ylabel('fp2', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
# ax2.set_ylim(5.0e-7, 1.0e-4)
ax2.legend(fontsize=20)

# plt.savefig("test.png")
plt.show()
plt.close(fig)