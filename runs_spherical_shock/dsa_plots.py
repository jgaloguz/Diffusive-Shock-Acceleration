# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import data
Nt_low = int(sys.argv[1])
Nt_high = int(sys.argv[2])
which = sys.argv[3]
t_arr = np.loadtxt("dsa_results/dsa_analytic_time.dat")
dsa_analytic_pos = np.loadtxt("dsa_results/dsa_analytic_pos_{:d}.dat".format(np.size(t_arr)-1))
dsa_analytic_mom = np.loadtxt("dsa_results/dsa_analytic_mom_{:d}.dat".format(np.size(t_arr)-1))
dsa_simulated_pos = [np.loadtxt("dsa_results/dsa_" + which + "_pos_{:d}_pp.dat".format(t)) for t in range(Nt_low, Nt_high)]
dsa_simulated_mom = [np.loadtxt("dsa_results/dsa_" + which + "_mom_{:d}_pp.dat".format(t)) for t in range(Nt_low, Nt_high)]

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax1 = fig.add_subplot(211, projection='rectilinear')

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
ax1.semilogy(dsa_analytic_pos[:,0], dsa_analytic_pos[:,1], color=colors[Nt_high-1])
for t in range(Nt_low, Nt_high):
   ax1.semilogy(dsa_simulated_pos[t-Nt_low][:,0], dsa_simulated_pos[t-Nt_low][:,1],
                color=colors[t], linestyle="", marker="o",
                label="t = {:.2e} day(s)".format(t_arr[t]))
ax1.set_xlabel('$x$ (au)', fontsize=20)
ax1.set_ylabel('$n$', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.legend(fontsize=20)

ax2 = fig.add_subplot(212, projection='rectilinear')

ax2.loglog(dsa_analytic_mom[:,0], dsa_analytic_mom[:,1], color=colors[Nt_high-1])
for t in range(Nt_low, Nt_high):
   ax2.loglog(dsa_simulated_mom[t-Nt_low][:,0], dsa_simulated_mom[t-Nt_low][:,1],
              color=colors[t], linestyle="", marker="s",
              label="t = {:.2e} day(s)".format(t_arr[t]))
ax2.set_xlabel('$E$ (MeV)', fontsize=20)
ax2.set_ylabel('$fp^2$', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.legend(fontsize=20)

# plt.savefig("test.png")
plt.show()
plt.close(fig)