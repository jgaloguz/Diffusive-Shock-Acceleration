# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import sys
from dsa_plot_common import InterpLikelihood, CondPathDens

# Load parameters
t = np.loadtxt("dsa_results/dsa_analytic_time.dat")
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")

# Load likelihood
hint = InterpLikelihood("dsa_results/feynman_kac/dsa_likelihood_space_time_alpha=1.00.dat", x)

# Import data from different trials
anlt0_x, anlt1_x, anlt2_x, anlt3_x, anlt0_n, anlt1_n, anlt2_n, anlt3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_analytic.dat", False, hint, x)
base0_x, base1_x, base2_x, base3_x, base0_n, base1_n, base2_n, base3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_base.dat", False, hint, x)
baseh0_x, baseh1_x, baseh2_x, baseh3_x, baseh0_n, baseh1_n, baseh2_n, baseh3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_base.dat", True, hint, x)
test0_x, test1_x, test2_x, test3_x, test0_n, test1_n, test2_n, test3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_ltest.dat", False, hint, x)
spltc0_x, spltc1_x, spltc2_x, spltc3_x, spltc0_n, spltc1_n, spltc2_n, spltc3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_0.dat", False, hint, x)
imps0_x, imps1_x, imps2_x, imps3_x, imps0_n, imps1_n, imps2_n, imps3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_imps_0.dat", False, hint, x)

# Report number of particles as a sanity check
print("Analytic solution number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(anlt0_n, anlt1_n, anlt2_n, anlt3_n))
print("Base method number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(base0_n, base1_n, base2_n, base3_n))
print("Splitting method number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(spltc0_n, spltc1_n, spltc2_n, spltc3_n))
print("Importance sampling number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(imps0_n, imps1_n, imps2_n, imps3_n))

# Plot solution
p_idx = 1
fig = plt.figure(figsize=(15, 10), layout='tight')

ax1 = fig.add_subplot(221, projection='rectilinear')

ax1.semilogy(x, anlt0_x, linewidth=3, label="analytic")
ax1.semilogy(x, base0_x, linewidth=3, linestyle="--", label="base method")
ax1.semilogy(x, baseh0_x, linewidth=3, label="base times likelihood")
ax1.semilogy(x, test0_x, linewidth=3, linestyle="--", label="base with reweighing")
ax1.semilogy(x, spltc0_x, linewidth=3, linestyle="-.", label="continuous splitting")
ax1.semilogy(x, imps0_x, linewidth=3, linestyle=":", label="importance sampling")
ax1.set_xlabel('$x$ (au)', fontsize=20)
ax1.set_ylabel('$n(x,t)$', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(-2.0,4.0)
ax1.set_ylim(1.0e-7,4.0e-3)
ax1.legend(fontsize=16)
ax1.annotate("$t=${:.0f}".format(t[0]), (3.0, 1.5e-7), fontsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

ax2.semilogy(x, anlt1_x, linewidth=3, label="analytic")
ax2.semilogy(x, base1_x, linewidth=3, linestyle="--", label="base method")
ax2.semilogy(x, baseh1_x, linewidth=3, label="base times likelihood")
ax2.semilogy(x, test1_x, linewidth=3, linestyle="--", label="base with reweighing")
ax2.semilogy(x, spltc1_x, linewidth=3, linestyle="-.", label="continuous splitting")
ax2.semilogy(x, imps1_x, linewidth=3, linestyle=":", label="importance sampling")
ax2.set_xlabel('$x$ (au)', fontsize=20)
ax2.set_ylabel('$n(x,t)$', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(-2.0,4.0)
ax2.set_ylim(1.0e-7,4.0e-3)
ax2.annotate("$t=${:.0f}".format(t[1]), (3.0, 1.5e-7), fontsize=20)

ax3 = fig.add_subplot(223, projection='rectilinear')

ax3.semilogy(x, anlt2_x, linewidth=3, label="analytic")
ax3.semilogy(x, base2_x, linewidth=3, linestyle="--", label="base method")
ax3.semilogy(x, baseh2_x, linewidth=3, label="base times likelihood")
ax3.semilogy(x, test2_x, linewidth=3, linestyle="--", label="base with reweighing")
ax3.semilogy(x, spltc2_x, linewidth=3, linestyle="-.", label="continuous splitting")
ax3.semilogy(x, imps2_x, linewidth=3, linestyle=":", label="importance sampling")
ax3.set_xlabel('$x$ (au)', fontsize=20)
ax3.set_ylabel('$n(x,t)$', fontsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(-2.0,4.0)
ax3.set_ylim(1.0e-7,4.0e-3)
ax3.annotate("$t=${:.0f}".format(t[2]), (3.0, 1.5e-7), fontsize=20)

ax4 = fig.add_subplot(224, projection='rectilinear')

ax4.semilogy(x, anlt3_x, linewidth=3, label="analytic")
ax4.semilogy(x, base3_x, linewidth=3, linestyle="--", label="base method")
ax4.semilogy(x, baseh3_x, linewidth=3, label="base time likelihood")
ax4.semilogy(x, test3_x, linewidth=3, linestyle="--", label="base with reweighing")
ax4.semilogy(x, spltc3_x, linewidth=3, linestyle="-.", label="continuous splitting")
ax4.semilogy(x, imps3_x, linewidth=3, linestyle=":", label="importance sampling")
ax4.set_xlabel('$x$ (au)', fontsize=20)
ax4.set_ylabel('$n(x,t)$', fontsize=20)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
ax4.set_xlim(-2.0,4.0)
ax4.set_ylim(1.0e-7,4.0e-3)
ax4.annotate("$t=${:.0f}".format(t[3]), (3.0, 1.5e-7), fontsize=20)

plt.savefig("dsa_results/dsa_splitting_continuous_likelihood_comparison.png", dpi=200)
plt.close(fig)