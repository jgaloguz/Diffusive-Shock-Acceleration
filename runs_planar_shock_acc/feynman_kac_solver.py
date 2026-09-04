# Calculate likelihood function using finite difference solver

# ============================================================
# Import libraries
# ============================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, eye
from scipy.sparse.linalg import splu
from scipy.interpolate import interp1d

# ============================================================
# Parameters (dimensionless)
# ============================================================

u1 = 4.0                # Upstream velocity
u2 = 1.0                # Downstream velocity
kappa0 = 0.668459       # Upstream diffusion coefficient
Lsh = 0.01              # Shock width
alpha = 1.0             # Splitting strength
tau = 1.73264           # Final time
tanh_fac = 4.0          # Thinning factor for the shock
Nt = 400000             # Number of time steps for solver
Nx = 10000              # Number of spatial grid points for solver
cluster = 0.02          # Cluster strength

# ============================================================
# Shock coefficients
# ============================================================

def stable_sech2(z):
    """
    Stable computation of sech(z)^2.

    Avoids overflow from np.cosh(z) when |z| is large.
    """
    a = np.exp(-2.0 * np.abs(z))
    return 4.0 * a / (1.0 + a)**2

def u_fun(x):
    return 0.5 * (u1 + u2) - 0.5 * (u1 - u2) * np.tanh(tanh_fac * x / Lsh)

def ux_fun(x):
    return -0.5 * (u1 - u2) * tanh_fac * stable_sech2(tanh_fac * x / Lsh) / Lsh

def kappa_fun(x):
    return kappa0 * u_fun(x)**2 / u1**2

def source_fun(x):
    """
    Source coefficient c(x) in

        g_s = L g + c g,

    where

        c(x) = -(alpha/3) u_x(x).
    """
    return -(alpha / 3.0) * ux_fun(x)

# ============================================================
# Shock-clustered nonuniform grid
# ============================================================

def make_shock_grid(x_left=-50.0, x_right=5.0, nx=1000, cluster_strength=0.5):
    """
    Construct a nonuniform grid clustered near x=0.

    The grid uses the map

        x = a sinh(s),

    with a comparable to Lsh.

    Smaller cluster_strength gives stronger clustering near the shock.
    """

    assert x_left < 0.0 < x_right
    assert nx >= 10

    a = cluster_strength * Lsh

    n_left = nx // 2 + 1
    n_right = nx - n_left + 1

    s_left = np.arcsinh(abs(x_left) / a)
    s_right = np.arcsinh(x_right / a)

    xL = -a * np.sinh(np.linspace(s_left, 0.0, n_left))
    xR =  a * np.sinh(np.linspace(0.0, s_right, n_right))

    # Drop duplicate zero from the right side
    x = np.concatenate([xL, xR[1:]])

    return x

# ============================================================
# Discrete backward-Kolmogorov operator
# ============================================================

def build_spatial_operator(x):
    """
    Build the spatial operator

        L f = u f_x + (kappa f_x)_x

    on a nonuniform grid.

    Important: this is the backward Kolmogorov operator. Therefore,
    for u > 0, the monotone upwind approximation is the forward
    difference

        u_i (f_{i+1} - f_i) / (x_{i+1} - x_i).

    This is opposite to the usual conservation-law upwind direction.
    """

    x = np.asarray(x)
    nx = len(x)

    u = u_fun(x)
    kappa = kappa_fun(x)

    L = lil_matrix((nx, nx))

    # Interior nodes
    for i in range(1, nx - 1):
        dxm = x[i] - x[i - 1]
        dxp = x[i + 1] - x[i]
        dxc = 0.5 * (dxm + dxp)

        # Conservative diffusion: (kappa f_x)_x
        km = 0.5 * (kappa[i - 1] + kappa[i])
        kp = 0.5 * (kappa[i] + kappa[i + 1])

        L[i, i - 1] += km / (dxc * dxm)
        L[i, i]     += -km / (dxc * dxm) - kp / (dxc * dxp)
        L[i, i + 1] += kp / (dxc * dxp)

        # Backward-Kolmogorov upwind advection: u f_x
        if u[i] >= 0.0:
            L[i, i + 1] += u[i] / dxp
            L[i, i]     += -u[i] / dxp
        else:
            L[i, i - 1] += -u[i] / dxm
            L[i, i]     += u[i] / dxm

    # Left boundary: reflecting/no-flux style closure
    dxp = x[1] - x[0]
    kp = 0.5 * (kappa[0] + kappa[1])

    L[0, 0] += -2.0 * kp / dxp**2
    L[0, 1] +=  2.0 * kp / dxp**2

    # Drift points into the domain only if u[0] > 0.
    # If it points outward, reflect by imposing zero derivative.
    if u[0] >= 0.0:
        L[0, 1] += u[0] / dxp
        L[0, 0] += -u[0] / dxp

    # Right boundary: reflecting/no-flux style closure
    dxm = x[-1] - x[-2]
    km = 0.5 * (kappa[-2] + kappa[-1])

    L[-1, -2] +=  2.0 * km / dxm**2
    L[-1, -1] += -2.0 * km / dxm**2

    # Drift points into the domain only if u[-1] < 0.
    # For your parameters u > 0, so this term is normally omitted
    # at the right boundary.
    if u[-1] <= 0.0:
        L[-1, -2] += -u[-1] / dxm
        L[-1, -1] += u[-1] / dxm

    return L.tocsr()

# ============================================================
# Feynman-Kac solver
# ============================================================

def solve_hbar(
    x_left=-50.0,
    x_right=10.0,
    nx=1000,
    nt=1000,
    cluster_strength=0.5,
    store_history=False,
):
    """
    Solve for hbar(t,x), especially hbar(0,x).

    We solve forward in s = tau - t:

        g_s = L g + c g,
        g(0,x) = 1,

    where

        c(x) = -(alpha/3) u_x(x).

    Then

        hbar(t,x) = g(tau - t, x).

    The source c is treated exactly by exponential multiplication,
    while L is treated by implicit Euler:

        g <- exp(dt c/2)
        g <- (I - dt L)^{-1} g
        g <- exp(dt c/2)

    This avoids the common NaN problem caused by treating the positive
    source term with the wrong implicit sign.
    """

    x = make_shock_grid(
        x_left=x_left,
        x_right=x_right,
        nx=nx,
        cluster_strength=cluster_strength,
    )

    dx = np.diff(x)
    print("Grid diagnostics:")
    print(f"  nx              = {nx}")
    print(f"  x_left, x_right = {x[0]:.6g}, {x[-1]:.6g}")
    print(f"  min dx          = {dx.min():.6e}")
    print(f"  max dx          = {dx.max():.6e}")
    print(f"  Lsh / min dx    = {Lsh / dx.min():.2f}")

    u = u_fun(x)
    ux = ux_fun(x)
    kappa = kappa_fun(x)
    c = source_fun(x)

    print("\nCoefficient diagnostics:")
    print(f"  min u, max u         = {u.min():.6g}, {u.max():.6g}")
    print(f"  min ux, max ux       = {ux.min():.6g}, {ux.max():.6g}")
    print(f"  min kappa, max kappa = {kappa.min():.6g}, {kappa.max():.6g}")
    print(f"  min c, max c         = {c.min():.6g}, {c.max():.6g}")

    if not np.all(np.isfinite(u)):
        raise FloatingPointError("u contains nonfinite values.")
    if not np.all(np.isfinite(ux)):
        raise FloatingPointError("ux contains nonfinite values.")
    if not np.all(np.isfinite(kappa)):
        raise FloatingPointError("kappa contains nonfinite values.")
    if not np.all(np.isfinite(c)):
        raise FloatingPointError("source c contains nonfinite values.")

    L = build_spatial_operator(x)
    dt = tau / nt

    # Stable implicit solve for the diffusion/advection generator L.
    # Since L has nonpositive diagonal and nonnegative off-diagonal
    # entries, I - dt L is well-conditioned in the M-matrix sense.
    M = eye(nx, format="csc") - dt * L
    lu = splu(M)

    half_source = np.exp(0.5 * dt * c)

    g = np.ones(nx)

    if store_history:
        G = np.zeros((nt // 10 + 1, nx))
        print("\nMemory usage:", G.size * 8 / 1.0e9, "GB")
        G[0, :] = g
    else:
        G = None

    for n in range(nt):
        g *= half_source
        g = lu.solve(g)
        g *= half_source

        if not np.all(np.isfinite(g)):
            raise FloatingPointError(
                f"Nonfinite values appeared at time step {n+1}. "
                f"Try increasing nt, reducing tau, or using a wider domain."
            )

        if store_history and (n+1)%10 == 0:
            G[n // 10 + 1, :] = g

    # g is hbar(0,x)
    hbar0 = g

    return x, hbar0, G

# ============================================================
# Run the solver
# ============================================================

print(f"alpha = {alpha: .2f}")
x, hbar0, G = solve_hbar(
    x_left=-50.0,
    x_right=10.0,
    nx=Nx,
    nt=Nt,
    cluster_strength=cluster,
    store_history=True,
)

# ============================================================
# Interpolator for hbar(0,x)
# ============================================================

hbar0_interp = interp1d(
    x,
    hbar0,
    kind="linear",
    bounds_error=False,
    fill_value=(hbar0[0], hbar0[-1]),
)

# Example: evaluate hbar(0,x) at selected positions
x_test = np.array([-40.0, -10.0, -1.0, -0.1, 0.0, 0.1, 1.0])
print("\nSample values of hbar(0,x):")
for xx, hh in zip(x_test, hbar0_interp(x_test)):
    print(f"  x = {xx: .4f}, hbar(0,x) = {hh:.8e}")

# ============================================================
# Plot and save
# ============================================================

print("\nPlotting and saving likelihood information...")

# Select times to plot and save likelihood
save_times = [1.0, 3.0, 10.0, 30.0]
save_idx =[int(Nt//10 * (1.0 - t / 30.0)) for t in save_times]
G_save = np.zeros((len(save_idx)+1, Nx))
G_save[0, :] = x[:]
for i in range(len(save_idx)):
    j = save_idx[i]
    G_save[i+1, :] = G[j, :]

plt.figure(figsize=(8, 4))
for i in range(len(save_idx)):
    plt.plot(G_save[0,:], G_save[i+1,:], label=f"t = {save_times[i]:.2f}")
plt.xlabel("x (au)")
plt.ylabel(r"$\overline{h}(t,x)$")
plt.grid(True)
plt.xlim(-12,2)
plt.tight_layout()
plt.legend()
plt.savefig("dsa_results/feynman_kac/dsa_likelihood_space_time_slices_alpha={:.2f}.png".format(alpha), dpi=300)

np.savetxt("dsa_results/feynman_kac/dsa_likelihood_space_time_alpha={:.2f}.dat".format(alpha), G_save)

# ============================================================
# Compute, plot, and output drift term 
# ============================================================

print("\nPlotting and saving induced drift information...")

# Initialize arrays
dtlogh = np.zeros((Nt // 10 + 1, Nx))
dxlogh = np.zeros((Nt // 10 + 1, Nx))
dx2logh = np.zeros((Nt // 10 + 1, Nx))
drift = np.zeros((Nt // 10 + 1, Nx))

# Iterate over time rows to compute derivatives
kappa = kappa_fun(x)
ux = ux_fun(x)
t = np.linspace(tau, 0.0, num=Nt // 10 + 1)
dt = t[1:] - t[:-1]
dx = x[1:] - x[:-1]
dx_pls = dx[1:]
dx_mns = dx[:-1]
dx2_pls = np.square(dx_pls)
dx2_mns = np.square(dx_mns)
denom = dx_pls * dx_mns * (dx_pls + dx_mns)
logh = np.log(G)
for i in range(Nt // 10 + 1):
    if i == 0:
        dtlogh[i,:] = (alpha / 3.0) * ux
    else:
        dtlogh[i,:] = (logh[i,:] - logh[i-1,:]) / dt[i-1]
    dxlogh[i,1:-1] = (dx2_mns * logh[i,2:] + (dx2_pls - dx2_mns) * logh[i,1:-1] - dx2_pls * logh[i,:-2]) / denom
    dx2logh[i,1:-1] = 2.0 * (dx_mns * logh[i,2:] - (dx_pls + dx_mns) * logh[i,1:-1] + dx_pls * logh[i,:-2]) / denom
    drift[i,1:-1] = 2.0 * kappa[1:-1] * dxlogh[i,1:-1]

# Select times to plot drift
drift_plot = np.zeros((len(save_idx), Nx))
for i in range(len(save_idx)):
    j = save_idx[i]
    drift_plot[i, :] = drift[j, :]

plt.figure(figsize=(8, 4))
for i in range(len(save_idx)):
    plt.plot(x, drift_plot[i,:], label=f"t = {save_times[i]:.2f}")
plt.plot(x, -2.0 * kappa * np.tanh(x / Lsh), label=f"Prinsloo")
plt.xlabel("x (au)")
plt.ylabel(r"2$\kappa(x)\partial_x\log\overline{h}(t,x)$ ($\times 10^7$ cm/s)")
plt.grid(True)
plt.xlim(-12,2)
plt.tight_layout()
plt.legend()
plt.savefig("dsa_results/feynman_kac/dsa_likelihood_drift_slices_alpha={:.2f}.png".format(alpha), dpi=300)

# Save only 1 out of every 10 time slices for efficiency
np.savetxt("dsa_results/feynman_kac/dsa_likelihood_x_alpha={:.2f}.dat".format(alpha), x)
np.savetxt("dsa_results/feynman_kac/dsa_likelihood_t_alpha={:.2f}.dat".format(alpha), t[::10])
np.savetxt("dsa_results/feynman_kac/dsa_likelihood_alpha={:.2f}.dat".format(alpha), G[::10,:])
np.savetxt("dsa_results/feynman_kac/dsa_likelihood_dt_alpha={:.2f}.dat".format(alpha), dtlogh[::10,:])
np.savetxt("dsa_results/feynman_kac/dsa_likelihood_dx_alpha={:.2f}.dat".format(alpha), dxlogh[::10,:])
np.savetxt("dsa_results/feynman_kac/dsa_likelihood_dx2_alpha={:.2f}.dat".format(alpha), dx2logh[::10,:])
