# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import sys
from dsa_plot_common import InterpLikelihood, CondPathDens

# Load parameters
t = np.loadtxt("dsa_results/dsa_analytic_time.dat")
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
conv = 1.496e+13 / 2.99792414e4 # conversion to physical units

# Load likelihood
hint = InterpLikelihood("dsa_results/feynman_kac/dsa_likelihood_space_time_alpha=1.00.dat", x)

# Import data from different trials
test0_x, test1_x, test2_x, test3_x, test0_n, test1_n, test2_n, test3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_ltest.dat", False, hint, x)
spltc0_x, spltc1_x, spltc2_x, spltc3_x, spltc0_n, spltc1_n, spltc2_n, spltc3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_split_0.dat", False, hint, x)
imps0_x, imps1_x, imps2_x, imps3_x, imps0_n, imps1_n, imps2_n, imps3_n \
   = CondPathDens("dsa_results/dsa_forward_path_dens_pp_imps_0.dat", False, hint, x)

# Report number of particles as a sanity check
print("Likelihood test number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(conv * test0_n, conv * test1_n, conv * test2_n, conv * test3_n))
print("Splitting method number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(conv * spltc0_n, conv * spltc1_n, conv * spltc2_n, conv * spltc3_n))
print("Importance sampling number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(conv * imps0_n, conv * imps1_n, conv * imps2_n, conv * imps3_n))

# Normalize
test1_x = test1_x / test1_n
spltc1_x = spltc1_x / spltc1_n
imps1_x = imps1_x / imps1_n
test3_x = test3_x / test3_n
spltc3_x = spltc3_x / spltc3_n
imps3_x = imps3_x / imps3_n

# Plot solution
threshold = 1.0e-3
fig = plt.figure(figsize=(15, 10), layout='tight')

ax1 = fig.add_subplot(221, projection='rectilinear')

ax1.semilogy(x, test1_x, linewidth=3, linestyle="-", label="terminal-weighted\ndirect paths")
ax1.semilogy(x, spltc1_x, linewidth=3, linestyle="--", label="continuous splitting")
ax1.semilogy(x, imps1_x, linewidth=3, linestyle=":", label="splitting-matched\nimportance sampling")
ax1.set_ylabel('$q_{t,m}(x)$ (au$^{-1}$)', fontsize=20)
ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(-2.0,4.0)
ax1.set_ylim(threshold,5.0e+0)
ax1.legend(fontsize=16)
ax1.set_title("$t=${:.0f} days".format(t[1]), fontsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

ax2.semilogy(x, test3_x, linewidth=3, linestyle="-", label="terminal-weighted\ndirect paths")
ax2.semilogy(x, spltc3_x, linewidth=3, linestyle="--", label="continuous splitting")
ax2.semilogy(x, imps3_x, linewidth=3, linestyle=":", label="splitting-matched\nimportance sampling")
ax2.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
ax2.set_xlim(-2.0,4.0)
ax2.set_ylim(threshold,5.0e+0)
ax2.set_title("$t=${:.0f} days".format(t[3]), fontsize=20)

# Mask arrays
test1_x[test1_x < threshold] = 0.0
spltc1_x[spltc1_x < threshold] = 0.0
imps1_x[imps1_x < threshold] = 0.0
test3_x[test3_x < threshold] = 0.0
spltc3_x[spltc3_x < threshold] = 0.0
imps3_x[imps3_x < threshold] = 0.0

ax3 = fig.add_subplot(223, projection='rectilinear')

ax3.semilogy(x, np.divide(np.abs(spltc1_x - test1_x), test1_x), linewidth=3, linestyle="--", label="continuous splitting", color="tab:orange")
ax3.semilogy(x, np.divide(np.abs(imps1_x - test1_x), test1_x), linewidth=3, linestyle=":", label="splitting-matched\nimportance sampling", color="tab:green")
ax3.set_xlabel('$x$ (au)', fontsize=20)
ax3.set_ylabel('$d_{t,m}(x)$', fontsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(-2.0,4.0)
ax3.set_ylim(1.0e-4,1.0e+0)
ax3.legend(fontsize=16)

ax4 = fig.add_subplot(224, projection='rectilinear')

ax4.semilogy(x, np.divide(np.abs(spltc3_x - test3_x), test3_x), linewidth=3, linestyle="--", label="continuous splitting", color="tab:orange")
ax4.semilogy(x, np.divide(np.abs(imps3_x - test3_x), test3_x), linewidth=3, linestyle=":", label="splitting-matched\nimportance sampling", color="tab:green")
ax4.set_xlabel('$x$ (au)', fontsize=20)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', which='both', left=False, labelleft=False)
ax4.set_xlim(-2.0,4.0)
ax4.set_ylim(1.0e-4,1.0e+0)

plt.savefig("dsa_results/dsa_splitting_continuous_likelihood_comparison_appendix.png", dpi=200)
plt.close(fig)