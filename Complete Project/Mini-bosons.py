import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from adjustText import adjust_text
from matplotlib.ticker import FuncFormatter
from PIL import Image
import io

# Constants (SI units)
G = 6.67430e-11
omega = 0.9  # Oscillation frequency (ω < 1 for bound state)
m = 1        # Scalar particle mass
r_inf=30
def solve_boson_star(omega, phi_center=0.01, r_center=1e-6, r_inf=30, num_points=500):
    def system(r, y, p):
        omega = p
        Phi, m_val, phi0, dphi0 = y
        eps = 1e-10
        r_adj = r + eps
        Lambda = -0.5 * np.log(np.clip(1 - 2 * m_val / r_adj, 1e-6, None))

        dLambda_dr = (4 * np.pi * r ** 3 * (
            omega ** 2 * phi0 ** 2 +
            np.exp(2 * (Phi - Lambda)) * dphi0 ** 2 +
            np.exp(2 * Phi) * phi0 ** 2
        ) - m_val) / (r * (r - 2 * m_val + eps))

        dm_dr = 4 * np.pi * r ** 2 * (
            omega ** 2 * phi0 ** 2 +
            np.exp(2 * (Phi - Lambda)) * dphi0 ** 2 +
            np.exp(2 * Phi) * phi0 ** 2
        )

        dPhi_dr = (
            (m_val + 4 * np.pi * r ** 3 * (
                dphi0 ** 2 +
                np.exp(2 * (Lambda - Phi)) * omega ** 2 * phi0 ** 2 -
                np.exp(2 * Lambda) * phi0 ** 2
            )) /
            (r * (r - 2 * m_val + eps))
        )

        ddphi0_dr = -(
            (2 / r + dLambda_dr - dPhi_dr) * dphi0 +
            (omega ** 2 * np.exp(2 * (Lambda - Phi)) - np.exp(2 * Lambda)) * phi0
        )

        return [dPhi_dr, dm_dr, dphi0, ddphi0_dr]

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

    r = np.linspace(r_center, r_inf, num_points)
    Phi_initial = np.zeros_like(r)
    m_initial = np.zeros_like(r)
    phi0_initial = phi_center * np.exp(-r ** 2)
    dphi0_initial = -2 * r * phi0_initial
    y_initial = np.vstack([Phi_initial, m_initial, phi0_initial, dphi0_initial])

    sol = solve_bvp(system, bc, r, y_initial, p=(omega,), max_nodes=10000)
    if not sol.success:
        print(f"Solver failed for ω = {omega}: {sol.message}")
        return None

    r_sol = sol.x
    Phi_sol, m_sol, phi0_sol, dphi0_sol = sol.y

    # Clip m(r) to non-negative values
    m_sol = np.clip(m_sol, 0, None)

    # Compute Lambda(r), avoiding division by zero
    Lambda_sol = -0.5 * np.log(np.clip(1 - 2 * m_sol / (r_sol + 1e-10), 1e-6, None))

    exp_neg_2Phi = np.exp(-2 * Phi_sol)
    exp_neg_2Lambda = np.exp(-2 * Lambda_sol)

    rho_sol = (omega ** 2 * exp_neg_2Phi * phi0_sol ** 2) + (exp_neg_2Lambda * dphi0_sol ** 2) + (m ** 2 * phi0_sol ** 2)
    P_r_sol = (exp_neg_2Lambda * dphi0_sol ** 2) + (omega ** 2 * exp_neg_2Phi * phi0_sol ** 2) - (m ** 2 * phi0_sol ** 2)
    P_perp_sol = (-exp_neg_2Lambda * dphi0_sol ** 2) + (m ** 2 * phi0_sol ** 2) - (omega ** 2 * exp_neg_2Phi * phi0_sol ** 2)

    return {
        'r': r_sol,
        'Phi': Phi_sol,
        'm': m_sol,
        'phi0': phi0_sol,
        'dphi0': dphi0_sol,
        'Lambda': Lambda_sol,
        'rho': rho_sol,
        'P_r': P_r_sol,
        'P_perp': P_perp_sol
    }

# Omega values
omega_values = np.round(np.linspace(0.1, 1, 50), 2)
m_omega = np.linspace(0.1, 1, 10)

frames_main = []
frames_pressure = []
frames_mass = []

def sci_notation(x, pos):
    if abs(x) < 1e-3 or abs(x) > 1e3:
        return f'{x:.2e}'
    else:
        return f'{x:.2f}'

r_inf_list = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

for r_inf in r_inf_list:
    frames_main = []
    frames_pressure = []
    masses = []

    for omega in omega_values:
        try:
            result = solve_boson_star(omega=omega, phi_center=0.01, r_center=1e-6, r_inf=r_inf, num_points=500)
            if result is None:
                continue
            '''
            # ---------- MAIN PROFILES ----------
            fig, axs = plt.subplots(2, 3, figsize=(18, 10))
            axs[0, 0].plot(result['r'], result['phi0'], 'k'); axs[0, 0].set_title(r'$\phi_0(r)$')
            axs[0, 1].plot(result['r'], result['Phi'], 'k'); axs[0, 1].set_title(r'$\Phi(r)$')
            axs[0, 2].plot(result['r'], result['Lambda'], 'k'); axs[0, 2].set_title(r'$\Lambda(r)$')

            mi = np.argmax(result['m'])
            r_clipped = result['r'][:mi + 1]
            m_clipped = np.clip(result['m'][:mi + 1], 0, None)
            axs[1, 0].plot(r_clipped, m_clipped, 'k', label=r'$m(r)$')
            axs[1, 0].set_title(f'$m(r)$ vs $r$')
            axs[1, 0].set_xlabel(r'$r$')
            axs[1, 0].set_ylabel(r'$m(r)$')

            axs[1, 1].plot(result['r'], result['rho'], 'k'); axs[1, 1].set_title(r'$\rho(r)$')
            axs[1, 2].plot(result['r'], result['P_r'], 'k', label=r'Radial Pressure $P_r$')
            axs[1, 2].plot(result['r'], result['P_perp'], 'k--', label=r'Tangential Pressure $P_\perp$')
            axs[1, 2].legend(); axs[1, 2].set_title(r'$P_r, P_\perp$')

            for ax in axs.flat:
                ax.set_xlabel('r')
                ax.grid(True)

            fig.suptitle(f'Boson Star Profiles — $\\omega = {omega:.3f}$, $r_{{\infty}} = {r_inf}$', fontsize=16)
            plt.tight_layout()
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=100)
            buf.seek(0)
            frames_main.append(Image.open(buf).copy())
            buf.close()
            plt.close(fig)
            '''
            # ---------- PRESSURE vs DENSITY ----------
            fig2, ax2 = plt.subplots(figsize=(8, 6))
            ax2.plot(result['rho'], result['P_r'], 'k', label=r'Radial Pressure $P_r$')
            ax2.plot(result['rho'], result['P_perp'], 'b--', label=r'Tangential Pressure $P_\perp$')
            ax2.set_xlabel(r'Energy Density $\rho$', fontsize=14)
            ax2.set_ylabel(r'Pressure', fontsize=14)
            ax2.set_title(rf'$P_r$ and $P_\perp$ vs $\rho$ for $\omega = {omega:.3f}$, $r_{{\infty}} = {r_inf}$', fontsize=16)
            ax2.xaxis.set_major_formatter(FuncFormatter(sci_notation))
            ax2.yaxis.set_major_formatter(FuncFormatter(sci_notation))
            ax2.legend(fontsize=12, loc='upper right')
            ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.tight_layout()

            buf2 = io.BytesIO()
            plt.savefig(buf2, format='png', dpi=100)
            buf2.seek(0)
            frames_pressure.append(Image.open(buf2).copy())
            buf2.close()
            plt.close(fig2)
            '''
            # ---------- RHO, Pr, P_perp vs RADIUS ----------
            r_plot = np.linspace(0, 10, len(result['rho']))

            fig3, ax3 = plt.subplots(figsize=(10, 6))
            ax3.plot(r_plot, result['rho'], 'k-', label=r'$\rho(r)$')
            ax3.plot(r_plot, result['P_r'], 'r--', label=r'$P_r(r)$')
            ax3.plot(r_plot, result['P_perp'], 'b:', label=r'$P_{\perp}(r)$')

            ax3.set_xlabel(r'$r$', fontsize=14)
            ax3.set_ylabel(r'Energy Density / Pressure', fontsize=14)
            ax3.set_title(fr'$\rho(r), P_r(r), P_\perp(r)$ vs r for $\omega = {omega:.2f}$, $r_\infty = {r_inf}$', fontsize=16)
            ax3.legend(fontsize=12)
            ax3.grid(True)
            plt.tight_layout()

            buf3 = io.BytesIO()
            plt.savefig(buf3, format='png', dpi=100)
            buf3.seek(0)
            frames_pressure.append(Image.open(buf3).copy())
            buf3.close()
            plt.close(fig3)
            '''
            print(f"Saved frame for ω = {omega:.2f}, r_inf = {r_inf}")

        except Exception as e:
            print(f"Error at ω = {omega:.2f}, r_inf = {r_inf}: {e}")

    # ---------- SAVE GIFs ----------
    if frames_main:
        frames_main[0].save(f"boson_star_profiles_rinf{r_inf}.gif",
                            save_all=True,
                            append_images=frames_main[1:],
                            duration=200,
                            loop=0,
                            optimize=False)
        print(f"GIF saved: 'boson_star_profiles_rinf{r_inf}.gif'")
    else:
        print(f"No main profile frames were generated for r_inf = {r_inf}")

    if frames_pressure:
        frames_pressure[0].save(f"pressure_vs_density_rinf{r_inf}.gif",
                                save_all=True,
                                append_images=frames_pressure[1:],
                                duration=200,
                                loop=0)
        print(f"GIF saved: 'pressure_vs_density_rinf{r_inf}.gif'")
    else:
        print(f"No pressure vs density frames were generated for r_inf = {r_inf}")
'''
for r_inf in r_inf_list:
    frames_main = []
    frames_pressure = []
    masses = []

    for omega in m_omega:
        try:
            result = solve_boson_star(omega=omega, phi_center=0.01, r_center=1e-6, r_inf=r_inf, num_points=500)
            if result is None:
                continue
            # ---------- ADM MASS vs OMEGA ----------
            if result is not None:
                M_ADM = result['m'][-1]
                masses.append((omega, M_ADM))
                if masses:
                    omegas_plot, masses_plot = zip(*masses)
                    plt.figure(figsize=(8, 6))
                    plt.plot(omegas_plot, masses_plot, 'k-', marker='o')
                    plt.xlabel(r'$\omega$', fontsize=14)
                    plt.ylabel(r'ADM Mass $M_{\text{ADM}}$', fontsize=14)
                    plt.title(fr'$M_{{\text{{ADM}}}}$ vs. $\omega$ for $r_{{\infty}} = {r_inf}$', fontsize=16)
                    plt.grid(True)

                    # Create labels
                    texts = []
                    for x, y in zip(omegas_plot, masses_plot):
                        texts.append(
                            plt.text(x, y, f'({x:.2f}, {y:.2f})', fontsize=9, ha='left', va='bottom')
                        )

                    # Run adjust_text with more iterations, and allow movement in x and y
                    adjust_text(
                        texts,
                        only_move={'points': 'y', 'text': 'xy'},
                        arrowprops=dict(arrowstyle='->', lw=0.5),
                        expand_text=(1.2, 1.4),
                        expand_points=(1.2, 1.4),
                        force_text=0.5,
                        force_points=0.2,
                        lim=200  # Increase iterations
                    )

                    plt.tight_layout()
                    plt.savefig(f"ADM_mass_vs_omega_rinf_{r_inf}.png")
                    plt.close()

            print(f"Saved frame for ω = {omega:.2f}, r_inf = {r_inf}")

        except Exception as e:
            print(f"Error at ω = {omega:.2f}, r_inf = {r_inf}: {e}")
'''