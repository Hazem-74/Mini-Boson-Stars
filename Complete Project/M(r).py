import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from matplotlib.ticker import FuncFormatter
from PIL import Image
from io import BytesIO

# Constants
G = 6.67430e-11
m = 1
phi_center = 0.01
r_center = 1e-6
r_inf_values = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
omega_values = np.linspace(0.1, 1, 10)

# Differential equations system
def system(r, y, p):
    omega = p
    Phi, m_func, phi0, dphi0 = y
    eps = 1e-10
    r_adj = r + eps

    Lambda = -0.5 * np.log(1 - 2 * m_func / r_adj)

    dLambda_dr = (4 * np.pi * r ** 3 * (
        omega ** 2 * phi0 ** 2 +
        np.exp(2 * (Phi - Lambda)) * dphi0 ** 2 +
        np.exp(2 * Phi) * phi0 ** 2
    ) - m_func) / (r * (r - 2 * m_func + eps))

    dm_dr = 4 * np.pi * r ** 2 * (
        omega ** 2 * phi0 ** 2 +
        np.exp(2 * (Phi - Lambda)) * dphi0 ** 2 +
        np.exp(2 * Phi) * phi0 ** 2
    )

    dPhi_dr = (
        (m_func + 4 * np.pi * r ** 3 * (
            dphi0 ** 2 +
            np.exp(2 * (Lambda - Phi)) * omega ** 2 * phi0 ** 2 -
            np.exp(2 * Lambda) * phi0 ** 2
        )) /
        (r * (r - 2 * m_func + eps))
    )

    ddphi0_dr = -(
        (2 / r + dLambda_dr - dPhi_dr) * dphi0 +
        (omega ** 2 * np.exp(2 * (Lambda - Phi)) - np.exp(2 * Lambda)) * phi0
    )

    return [dPhi_dr, dm_dr, dphi0, ddphi0_dr]

# Boundary conditions
def bc(ya, yb, p):
    Phi_a, m_a, phi0_a, dphi0_a = ya
    Phi_b, m_b, phi0_b, dphi0_b = yb
    return [
        phi0_a - phi_center,
        dphi0_a,
        phi0_b,
        Phi_b,
        m_a
    ]

# Plotting function (in memory)
def plot_mass_profile_to_image(r_vals, m_vals, omega, r_inf):
    mi = np.argmax(m_vals)
    max_r = r_vals[mi]

    plt.figure(figsize=(8, 6))
    plt.plot(r_vals, np.clip(m_vals, 0, None), 'k-', label=r'$m(r)$')
    plt.title(f'Mass Profile $m(r)$ vs Radius $r$\n$\omega={omega:.3f}$, $r_{{\infty}}={r_inf}$')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$m(r)$')
    plt.xlim(0, max_r)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    buf = BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    buf.seek(0)
    return Image.open(buf).convert('P', palette=Image.ADAPTIVE)

# Outer loop over different r_inf values
for r_inf in r_inf_values:
    print(f"\n=== Solving for r_inf = {r_inf} ===")

    # Redefine r and initial guess for this r_inf
    r = np.linspace(r_center, r_inf, 1000)
    Phi_initial = np.zeros_like(r)
    m_initial = np.zeros_like(r)
    phi0_initial = phi_center * np.exp(-r ** 2)
    dphi0_initial = -2 * r * phi_center * np.exp(-r ** 2)
    y_initial = np.vstack([Phi_initial, m_initial, phi0_initial, dphi0_initial])

    frames = []

    print("Generating frames...")
    for i, omega in enumerate(omega_values):
        print(f"  Solving for omega = {omega:.3f}")
        sol = solve_bvp(system, bc, r, y_initial, p=(omega,), tol=1e-3, max_nodes=1000)

        if sol.status != 0:
            print(f"  Warning: solver did not converge for omega = {omega:.3f}")

        r_sol = sol.x
        m_sol = sol.y[1]
        frame = plot_mass_profile_to_image(r_sol, m_sol, omega, r_inf)
        frames.append(frame)

    # Create GIF for this r_inf
    print(f"Creating GIF for r_inf = {r_inf}...")
    gif_filename = f"mass_profile_vs_omega_rinf{r_inf}.gif"
    frames[0].save(gif_filename,
                   save_all=True,
                   append_images=frames[1:],
                   duration=300,
                   loop=0)
    print(f"GIF saved as '{gif_filename}'")