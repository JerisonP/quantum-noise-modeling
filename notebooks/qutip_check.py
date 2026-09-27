"""
qutip_check.py: solve the thesis model with QuTiP (an independent, widely used
quantum-dynamics package) on the SAME noise paths the Julia package used, and
compare trajectory by trajectory.

Needs:  pip install qutip numpy
Input (written by export_for_qutip.jl):
  results/qutip_noise.csv    the noise paths   (package CSV format: t, xi_1 … xi_M)
  results/package_out.csv    the package's per-trajectory ε and final states
"""
import sys, time
import numpy as np
import qutip as qt

def read_table(path):
    """Minimal reader for the package's CSV format: skip '#' lines, one header row."""
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    cols = lines[0].split(",")
    data = np.array([[float(x) for x in l.split(",")] for l in lines[1:]])
    return cols, data

theta, tg = np.pi / 2, 1.0                                   # thesis: θ = π/2, t_g = τ_c = 1
cols, noise = read_table("results/qutip_noise.csv")
t = noise[:, cols.index("t")]
xis = noise[:, 1:]                                           # one column per trajectory
_, pkg = read_table("results/package_out.csv")               # columns: eps, zp_rho11, zp_rho22, zp_re_rho12, zp_im_rho12, xp_...

# Thesis Eq. 3.12, built from QuTiP's own operators: H(t) = f_x(t)/2 σx + ξ(t)/2 σz
f = lambda s: theta / tg * (1 - np.cos(2 * np.pi * s / tg))
Uideal = (-1j * theta / 2 * qt.sigmax()).expm()              # Eq. 3.8
opts = {"atol": 1e-12, "rtol": 1e-10, "max_step": t[1] - t[0], "nsteps": 10**6}

M = xis.shape[1] if len(sys.argv) < 2 else int(sys.argv[1])
eps, rho_z, rho_x = [], [], []
psi_z, psi_x = qt.basis(2, 0), (qt.basis(2, 0) + qt.basis(2, 1)).unit()
t0 = time.time()
for j in range(M):
    xi = lambda s, j=j: np.interp(s, t, xis[:, j])           # the same linear interpolation
    H = qt.QobjEvo([[0.5 * qt.sigmax(), f], [0.5 * qt.sigmaz(), xi]])
    U = qt.sesolve(H, qt.qeye(2), [0, tg], options=opts).states[-1]   # propagator U(t_g)
    eps.append(1 - (abs((Uideal.dag() * U).tr()) ** 2 + 2) / 6)       # Nielsen, d = 2
    rho_z.append((U * psi_z).proj().full())
    rho_x.append((U * psi_x).proj().full())
eps = np.array(eps)
print(f"QuTiP {qt.__version__}: {M} trajectories in {time.time() - t0:.1f} s")
print(f"  eps (QuTiP)   = {eps.mean():.6f} ± {eps.std(ddof=1) / np.sqrt(M):.6f}")
print(f"  eps (package) = {pkg[:M, 0].mean():.6f}")
print(f"  largest per-trajectory |eps_QuTiP - eps_package| = {np.abs(eps - pkg[:M, 0]).max():.1e}   (pass: < 1e-6)")
for name, rhos, c in (("|0>", rho_z, 1), ("|+>", rho_x, 5)):
    q = np.mean(rhos, axis=0)
    p = pkg[:M, c:c + 4].mean(axis=0)                        # rho11, rho22, Re rho12, Im rho12
    d = max(abs(q[0, 0].real - p[0]), abs(q[1, 1].real - p[1]), abs(q[0, 1].real - p[2]), abs(q[0, 1].imag - p[3]))
    print(f"  final <rho> from {name}: largest element difference = {d:.1e}   (pass: < 1e-6)")