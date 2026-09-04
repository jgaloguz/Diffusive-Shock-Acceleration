# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import sys
from dsa_plot_common import InterpLikelihood, CondPathDens

# Load parameters
alpha = 1.0
t = np.loadtxt("dsa_results/dsa_analytic_time.dat")
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")

hint = InterpLikelihood("dsa_results/feynman_kac/dsa_likelihood_space_time_alpha=1.00.dat", x)

spltc0_x, spltc1_x, spltc2_x, spltc3_x, spltc0_n, spltc1_n, spltc2_n, spltc3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_0.dat", False, hint, x)
spltd10_x, spltd11_x, spltd12_x, spltd13_x, spltd10_n, spltd11_n, spltd12_n, spltd13_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_10.dat", False, hint, x)
spltd20_x, spltd21_x, spltd22_x, spltd23_x, spltd20_n, spltd21_n, spltd22_n, spltd23_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_100.dat", False, hint, x)
spltd30_x, spltd31_x, spltd32_x, spltd33_x, spltd30_n, spltd31_n, spltd32_n, spltd33_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_1000.dat", False, hint, x)
spltd40_x, spltd41_x, spltd42_x, spltd43_x, spltd40_n, spltd41_n, spltd42_n, spltd43_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_10000.dat", False, hint, x)

print("Continuous splitting method number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(spltc0_n, spltc1_n, spltc2_n, spltc3_n))
print("Discrete splitting method (10 thresholds) number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(spltd10_n, spltd11_n, spltd12_n, spltd13_n))
print("Discrete splitting method (100 thresholds) number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(spltd20_n, spltd21_n, spltd22_n, spltd23_n))
print("Discrete splitting method (1000 thresholds) number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(spltd30_n, spltd31_n, spltd32_n, spltd33_n))
print("Discrete splitting method (10000 thresholds) number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(spltd40_n, spltd41_n, spltd42_n, spltd43_n))

# Plot solution
p_idx = 1
fig = plt.figure(figsize=(15, 10), layout='tight')

ax1 = fig.add_subplot(221, projection='rectilinear')

ax1.semilogy(x, spltc0_x, linewidth=3, label="continuous splitting")
ax1.semilogy(x, spltd10_x, linewidth=3, linestyle="-", label="discrete splitting (10)")
ax1.semilogy(x, spltd20_x, linewidth=3, linestyle="--", label="discrete splitting (100)")
ax1.semilogy(x, spltd30_x, linewidth=3, linestyle="-.", label="discrete splitting (1000)")
ax1.semilogy(x, spltd40_x, linewidth=3, linestyle=":", label="discrete splitting (10000)")
ax1.set_xlabel('$x$ (au)', fontsize=20)
ax1.set_ylabel('$n(x,t)$', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(-2.0,4.0)
ax1.set_ylim(1.0e-7,4.0e-3)
ax1.legend(fontsize=16)
ax1.annotate("$t=${:.0f}".format(t[0]), (3.0, 1.5e-7), fontsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

ax2.semilogy(x, spltc1_x, linewidth=3, label="continuous splitting")
ax2.semilogy(x, spltd11_x, linewidth=3, linestyle="-", label="discrete splitting (10)")
ax2.semilogy(x, spltd21_x, linewidth=3, linestyle="--", label="discrete splitting (100)")
ax2.semilogy(x, spltd31_x, linewidth=3, linestyle="-.", label="discrete splitting (1000)")
ax2.semilogy(x, spltd41_x, linewidth=3, linestyle=":", label="discrete splitting (10000)")
ax2.set_xlabel('$x$ (au)', fontsize=20)
ax2.set_ylabel('$n(x,t)$', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(-2.0,4.0)
ax2.set_ylim(1.0e-7,4.0e-3)
ax2.annotate("$t=${:.0f}".format(t[1]), (3.0, 1.5e-7), fontsize=20)

ax3 = fig.add_subplot(223, projection='rectilinear')

ax3.semilogy(x, spltc2_x, linewidth=3, label="continuous splitting")
ax3.semilogy(x, spltd12_x, linewidth=3, linestyle="-", label="discrete splitting (10)")
ax3.semilogy(x, spltd22_x, linewidth=3, linestyle="--", label="discrete splitting (100)")
ax3.semilogy(x, spltd32_x, linewidth=3, linestyle="-.", label="discrete splitting (1000)")
ax3.semilogy(x, spltd42_x, linewidth=3, linestyle=":", label="discrete splitting (10000)")
ax3.set_xlabel('$x$ (au)', fontsize=20)
ax3.set_ylabel('$n(x,t)$', fontsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(-2.0,4.0)
ax3.set_ylim(1.0e-7,4.0e-3)
ax3.annotate("$t=${:.0f}".format(t[2]), (3.0, 1.5e-7), fontsize=20)

ax4 = fig.add_subplot(224, projection='rectilinear')

ax4.semilogy(x, spltc3_x, linewidth=3,label="continuous splitting")
ax4.semilogy(x, spltd13_x, linewidth=3, linestyle="-", label="discrete splitting (10)")
ax4.semilogy(x, spltd23_x, linewidth=3, linestyle="--", label="discrete splitting (100)")
ax4.semilogy(x, spltd33_x, linewidth=3, linestyle="-.", label="discrete splitting (1000)")
ax4.semilogy(x, spltd43_x, linewidth=3, linestyle=":", label="discrete splitting (10000)")
ax4.set_xlabel('$x$ (au)', fontsize=20)
ax4.set_ylabel('$n(x,t)$', fontsize=20)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
ax4.set_xlim(-2.0,4.0)
ax4.set_ylim(1.0e-7,4.0e-3)
ax4.annotate("$t=${:.0f}".format(t[3]), (3.0, 1.5e-7), fontsize=20)

plt.savefig("dsa_results/dsa_splitting_continuous_discrete_comparison.png", dpi=200)
plt.close(fig)