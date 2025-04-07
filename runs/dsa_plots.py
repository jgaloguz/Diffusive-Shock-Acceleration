# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Import data
Nt = 5
t_arr = np.loadtxt("dsa_results/dsa_analytic_time.dat")
dsa_analytic_pos = [np.loadtxt("dsa_results/dsa_analytic_pos_{:d}.dat".format(t)) for t in range(Nt)]
dsa_analytic_mom = [np.loadtxt("dsa_results/dsa_analytic_mom_{:d}.dat".format(t)) for t in range(Nt)]
dsa_forward_pos = [np.loadtxt("dsa_results/dsa_forward_pos_{:d}_pp.dat".format(t)) for t in range(Nt)]
dsa_forward_mom = [np.loadtxt("dsa_results/dsa_forward_mom_{:d}_pp.dat".format(t)) for t in range(Nt)]
dsa_backward_mom = [np.loadtxt("dsa_results/dsa_backward_mom_{:d}_pp.dat".format(t)) for t in range(Nt)]

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax1 = fig.add_subplot(211, projection='rectilinear')

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
for t in range(Nt):
   ax1.semilogy(dsa_analytic_pos[t][:,0], dsa_analytic_pos[t][:,1],
                color=colors[t], label="t = {:.2e}".format(t_arr[t]))
   ax1.semilogy(dsa_forward_pos[t][:,0], dsa_forward_pos[t][:,1],
                color=colors[t], linestyle="", marker="o")
ax1.set_xlabel('x', fontsize=20)
ax1.set_ylabel('N', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(-1.0, 4.0)
ax1.set_ylim(1.0e-4, 2.0e0)
ax1.legend(fontsize=20)

ax2 = fig.add_subplot(212, projection='rectilinear')

for t in range(Nt):
   ax2.loglog(dsa_analytic_mom[t][:,0], dsa_analytic_mom[t][:,1],
              color=colors[t], label="t = {:.2e}".format(t_arr[t]))
   ax2.loglog(dsa_backward_mom[t][:,0], dsa_backward_mom[t][:,1],
              color=colors[t], linestyle="", marker="s")
ax2.set_xlabel('p', fontsize=20)
ax2.set_ylabel('fp2', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_ylim(1.0e-8, 2.0e-4)
ax2.legend(fontsize=20)

# plt.savefig("test.png")
plt.show()
plt.close(fig)