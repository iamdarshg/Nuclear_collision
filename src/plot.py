from plot import *
# import plot
import os
import matplotlib.pyplot as plt
import numpy as np
import platform

# --- Configuration ---
# Determine the correct library filename based on the OS
# if platform.system() == "Linux" or platform.system() == "Darwin": # Linux or macOS
#     LIB_FILENAME = "libeventsimulator.so"
# else: # Assume Windows
#     LIB_FILENAME = "libeventsimulator.a"
# os.add_dll_directory(os.path.dirname(__file__))  # Ensure the directory is added to the DLL search path
# # Path to the shared library
# # Ensure this script is in the same directory as the compiled C++ library
# LIB_PATH = os.path.join(os.path.dirname(__file__), LIB_FILENAME)

# # Simulation parameters for Python side
NUM_ENERGY_POINTS = 200           # Number of different energy values to test



# # --- Load C++ Shared Library ---
# try:
#     # Load the shared library
#     event_simulator_lib = ctypes.CDLL(LIB_PATH)

#     # Define the argument types and return type of the C++ function
#     # bool simulate_collision_event(double v_initial, double impact_parameter_val)
#     event_simulator_lib.simulate_collision_event.argtypes = [ctypes.c_double, ctypes.c_double]
#     event_simulator_lib.simulate_collision_event.restype = ctypes.c_bool

#     print(f"Successfully loaded C++ library: {LIB_PATH}")

# except OSError as e:
#     print(f"Error loading C++ library: {e}")
#     print(f"Please ensure '{LIB_FILENAME}' is compiled and located in the same directory as this script.")
#     print("Compilation commands:")
#     print("  Linux/macOS: g++ -shared -o libeventsimulator.so -fPIC event_simulator.cpp")
#     print("  Windows:     g++ -shared -o libeventsimulator.dll event_simulator.cpp")
#     exit()

# --- Simulation Logic in Python ---
event_simulator_lib = None
def calculate_collision_probability(energy_vel):
    """
    Calculates the collision probability for a given energy.
    This involves running multiple Monte Carlo trials over a range of impact parameters.
    """
    is_collision = calculate_collision_probabigy_MeV(energy_vel)
    return is_collision

# --- Main Plotting Execution ---

if __name__ == "__main__":

    # min_energy_MeV = 1 # MeV
    # max_energy_MeV = 1000 # MeV
    # print(dir())
    # energies_MeV = np.linspace(min_energy_MeV, max_energy_MeV, NUM_ENERGY_POINTS)
    # energies_joules = energies_MeV * 1.602e-13 # Convert MeV to Joules

    # collision_probabilities = []

    # print("\nStarting simulation for different energies...")
    # for i, energy_j in enumerate(energies_joules):
    #     print(f"  Simulating for Energy: {energies_MeV[i]:.2f} MeV ({energy_j:.2e} J)...")
    #     prob = calculate_collision_probability(energy_j)
    #     collision_probabilities.append(prob)
    #     print(f"    Collision Probability: {prob:.4f}")

    # print("\nSimulation complete. Plotting results...")

    # # --- Plotting ---
    # plt.figure(figsize=(10, 6))
    # plt.plot(energies_MeV, collision_probabilities, marker='o', linestyle='-', color='teal', label='Collision Probability')

    # plt.title(f'{os.path.basename(LIB_FILENAME).replace(".so", "").replace(".dll", "")}: Collision Probability vs. Energy', fontsize=16)
    # plt.xlabel('Total Kinetic Energy (MeV)', fontsize=14)
    # plt.ylabel('Collision Probability', fontsize=14)
    # plt.grid(True, linestyle='--', alpha=0.7)
    # plt.ylim(-0.05, 1.05) # Probability is between 0 and 1
    # plt.legend(fontsize=12)
    # plt.tick_params(axis='both', which='major', labelsize=12)
    # plt.tight_layout()

    # # Save the plot
    # plot_filename = "collision_probability_vs_energy.png"
    # plt.savefig(plot_filename, dpi=300)
    # print(f"\nPlot saved as '{plot_filename}'")

    # # Show the plot
    # plt.show()
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import simpson
    from scipy.special import spherical_jn

    N = 120
    X, Y, Z = np.mgrid[-25:25:N*1j, -25:25:N*1j, -25:25:N*1j]
    r = np.sqrt(X**2 + Y**2 + Z**2)

    # Define the 1D radial wavefunctions
    def u_s_wave(r_fm):
        return r_fm * np.exp(-0.7 * r_fm) * (1 + 0.2 * r_fm)

    def w_d_wave(r_fm):
        return r_fm**3 * np.exp(-1.6 * r_fm) * (1 + 0.3 * r_fm)

    # Convert the 3D grid 'r' (in a0 units) to fm for the wavefunction definitions
    # This conversion factor is approximate, a more precise one might be needed
    a0_in_fm = 0.529 * 1e-5 / 1e-15
    r_fm = r * a0_in_fm

    # Create the initial 3D wavefunctions
    psi_S_initial = u_s_wave(r_fm)
    psi_D_initial = w_d_wave(r_fm)

    # Normalize them separately
    psi_S_initial /= np.sqrt(np.sum(np.abs(psi_S_initial)**2))
    psi_D_initial /= np.sqrt(np.sum(np.abs(psi_D_initial)**2))

    def simulate_deuteron_potential(waveform, plot=False):
        """
        Simulates the potential well of a deuteron given a stacked waveform.

        Args:
            waveform (function): A function that takes r (in fm) and returns a tuple
                                of S- and D-wave radial wavefunctions.

        Returns:
            tuple: A tuple containing the expectation values for the tensor, central, and total potentials.
        """
        # Constants
        f_pi=92.4; g_A=1.29; m_pi=138.0; hbar_c=197.3
        Lambda=500.0; Lambda_fm = Lambda/hbar_c
        q_max=1000

        r_vals_fm = np.linspace(0.01, 20.0, 1000)
        q_vals = np.linspace(1e-4, 2000, 1000)

        # Deuteron: S=1, T=0
        tau_dot = -3  # <τ₁·τ₂> = -3 for T=0
        sigma_dot = 1 # <σ₁·σ₂> = +1
        # Tensor operator average for S-D mixing: use approximate <S12> ~ -2 for deuteron major component
        S12_avg = -2

        def V_OPE(r):
            pref = (g_A**2)/(4*np.pi*f_pi**2)
            W_S = m_pi**2 * np.exp(-m_pi*r)/r
            W_T = (1 + 3/(m_pi*r) + 3/(m_pi*r)**2) * np.exp(-m_pi*r)/r
            return tau_dot * pref * (sigma_dot * W_S + S12_avg * W_T)

        def V_contact(r):
            delta = (Lambda_fm**3/np.pi**1.5) * np.exp(-(Lambda_fm*r)**2)
            C_S, C_T = 400, -100
            return C_S * delta + C_T * delta * sigma_dot

        def V_TPE_momentum(q_mu):
            return q_mu**2 * np.exp(-q_mu**2 / Lambda**2)

        def V_TPE_bessel(r, q_vals):
            integrand = q_vals**2 * V_TPE_momentum(np.sqrt(q_vals**2 + m_pi**2)) * spherical_jn(0, np.outer(r, q_vals))
            integral = simpson(integrand, q_vals, axis=-1)
            return (2/np.pi) * integral

        def V_tensor(r):
            # simplified spatial profile of tensor force strength
            return -50 * np.exp(-r/1.0) / r**3  # arbitrary shape

        def V_total(r, q_vals):
            return V_OPE(r) + V_contact(r) + V_TPE_bessel(r, q_vals) + S12_avg * V_tensor(r)

        # Normalize the wavefunctions
        def normalize_wf(u, w, r):
            norm = simpson(u**2 + w**2, r)
            return u / np.sqrt(norm), w / np.sqrt(norm)

        u_raw, w_raw = waveform(r_vals_fm)
        u_norm, w_norm = normalize_wf(u_raw, w_raw, r_vals_fm)

        # Tensor operator S12 expectation (approximation)
        # <S12> for deuteron averaged over spin orientation is ~0.5 for S-D mixing
        S12 = 0.5


        g_A = 1.29
        f_pi = 92.4  # MeV
        m_pi = 138.0  # MeV
        hbar_c = 197.3  # MeV·fm

        # Loop function L(q) in TPE
        def loop_L(q):
            sqrt_term = np.sqrt(4 * m_pi**2 + q**2)
            return (sqrt_term / q) * np.log((sqrt_term + q) / (sqrt_term - q) + 1e-12)

        # V_C_TPE in momentum space (NLO)
        def V_C_TPE_q(q):
            factor = -3 * g_A**4 / (64 * np.pi**2 * f_pi**4)
            return factor * loop_L(q) * (4 * m_pi**2 + q**2)

        # Coordinate-space central TPE potential
        def V_C_TPE_r(r, q_vals):
            integrand = q_vals**2 * V_C_TPE_q(q_vals) * spherical_jn(0, np.outer(r, q_vals))
            integral = simpson(integrand, q_vals, axis=-1)
            return integral / (2 * np.pi**2)


        # Compute expectation values using S- and D-wave functions
        def expectation_value(V, u, w, r):
            return simpson(u**2 * V + w**2 * V + 2 * u * w * S12 * V, r)

        V_ten_vals = V_tensor(r_vals_fm)
        V_cen_vals = V_C_TPE_r(r_vals_fm, q_vals)

        E_tensor = expectation_value(V_ten_vals, u_norm, w_norm, r_vals_fm)
        E_central = expectation_value(V_cen_vals, u_norm, w_norm, r_vals_fm)
        E_total = E_tensor + E_central
        if plot:
        # Plot the wavefunctions and energy densities
            plt.figure(figsize=(10, 5))
            plt.subplot(1, 2, 1)
            plt.plot(r_vals_fm, u_norm, label='S-wave (u)')
            plt.plot(r_vals_fm, w_norm, label='D-wave (w)')
            plt.title('Deuteron Wavefunctions')
            plt.xlabel('r (fm)')
            plt.ylabel('u(r), w(r)')
            plt.grid(True)
            plt.legend()

            plt.subplot(1, 2, 2)
            plt.plot(r_vals_fm, u_norm**2 * V_ten_vals, label='Tensor Energy Density')
            # plt.plot(r_vals_fm, u_norm**2 * V_cen_vals, label='Central Energy Density')
            plt.plot(r_vals_fm, u_norm**2 * (V_cen_vals+V_cen_vals), label='Total Energy Density')
            plt.title('Deuteron Energy Densities')
            plt.xlabel('r (fm)')
            plt.ylabel('Energy Density (MeV/fm)')
            plt.yscale('log')
            plt.grid(True)
            plt.legend()

            plt.tight_layout()
            plt.show()

        return (E_tensor, E_central, E_total)

    # Example usage with the original wavefunctions
    def stacked_waveform(r):
        u_s = r * np.exp(-0.7 * r) * (1 + 0.2 * r)
        w_d = r**3 * np.exp(-1.6 * r) * (1 + 0.3 * r)
        return u_s, w_d

    # simulate_deuteron_potential(stacked_waveform)


    def get_deuteron_potential(waveform, r_vals_fm):
        return simulate_deuteron_potential(waveform)

    def get_central_potential(r_fm):
        """Returns a simplified central potential."""
        # A simple attractive Gaussian well for demonstration
        V0 = -50  # MeV
        a = 2.0   # fm
        return V0 * np.exp(-(r_fm**2) / a**2)

    def get_tensor_potential(r_fm):
        """Returns a simplified tensor potential for coupling."""
        # This is a schematic form. A realistic tensor potential is more complex.
        V_tensor_strength = -30 # MeV
        return V_tensor_strength * np.exp(-r_fm / 1.5)

    # Create 3D potential fields from the 1D definitions
    V_central = get_central_potential(r_fm)
    V_tensor = get_tensor_potential(r_fm)

    # Step 3 & 4: Create and Implement the Coupled Evolution Function

    def split_step_coupled(psi_S, psi_D, V_central, V_tensor, dt, mask):
        """
        Evolves two coupled wavefunctions (S and D) over one time step.

        Args:
            psi_S (np.ndarray): The S-wave component.
            psi_D (np.ndarray): The D-wave component.
            V_central (np.ndarray): The central potential.
            V_tensor (np.ndarray): The tensor (coupling) potential.
            dt (float): The time step.
            mask (np.ndarray): The absorbing boundary mask.

        Returns:
            tuple: A tuple containing the updated psi_S and psi_D.
        """
        # --- Potential Step (Part 1) ---
        # Apply central potential to both
        psi_S_half = np.exp(-0.5j * V_central * dt) * psi_S
        psi_D_half = np.exp(-0.5j * V_central * dt) * psi_D

        # Apply coupling potential
        # This is a simplified coupling scheme. A more rigorous one would use matrix exponentiation.
        psi_S_coupled = psi_S_half - 0.5j * V_tensor * dt * psi_D_half
        psi_D_coupled = psi_D_half - 0.5j * V_tensor * dt * psi_S_half


        # --- Kinetic Step ---
        # Transform to momentum space
        psi_S_k = np.fft.fftn(psi_S_coupled)
        psi_D_k = np.fft.fftn(psi_D_coupled)

        # Define momentum grid
        kx = np.fft.fftfreq(psi_S.shape[0])
        ky = np.fft.fftfreq(psi_S.shape[1])
        kz = np.fft.fftfreq(psi_S.shape[2])
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        K2 = KX**2 + KY**2 + KZ**2 # Kinetic energy term

        # Evolve in momentum space
        psi_S_k_evolved = np.exp(-1j * K2 * dt) * psi_S_k
        psi_D_k_evolved = np.exp(-1j * K2 * dt) * psi_D_k

        # Transform back to real space
        psi_S_evolved = np.fft.ifftn(psi_S_k_evolved)
        psi_D_evolved = np.fft.ifftn(psi_D_k_evolved)


        # --- Potential Step (Part 2) ---
        # Apply central potential again
        psi_S_final_half = np.exp(-0.5j * V_central * dt) * psi_S_evolved
        psi_D_final_half = np.exp(-0.5j * V_central * dt) * psi_D_evolved

        # Apply coupling potential again
        psi_S_final = psi_S_final_half - 0.5j * V_tensor * dt * psi_D_final_half
        psi_D_final = psi_D_final_half - 0.5j * V_tensor * dt * psi_S_final_half

        # --- Apply Absorbing Mask ---
        psi_S_final *= mask
        psi_D_final *= mask

        return psi_S_final, psi_D_final

    # Step 5: Update the Main Loop and Plot
    def get_absorbing_boundary_mask(shape, width, strength):
        """
        Creates an absorbing boundary mask.

        Args:
            shape (tuple): The shape of the mask.
            width (int): The width of the absorbing boundary.
            strength (float): The strength of the absorbing boundary.

        Returns:
            np.ndarray: The absorbing boundary mask.
        """
        mask = np.ones(shape)
        for i in range(3):
            for j in range(width):
                mask[j,:,:] *= (1 - strength * (1 - j/width)**2)
                mask[-j-1,:,:] *= (1 - strength * (1 - j/width)**2)
                mask[:,j,:] *= (1 - strength * (1 - j/width)**2)
                mask[:,-j-1,:] *= (1 - strength * (1 - j/width)**2)
                mask[:,:,j] *= (1 - strength * (1 - j/width)**2)
                mask[:,:,-j-1] *= (1 - strength * (1 - j/width)**2)
        return mask
    from plotly import graph_objects as go
    # Time evolution parameters
    dt = 0.00001
    num_steps = 3000
    plot_every = 20

    # Get the potential using the function we defined earlier
    r_vals_fm = np.linspace(0.01, 20.0, N)
    def waveform(r):
        u_s = r * np.exp(-0.7 * r) * (1 + 0.2 * r)
        w_d = r**3 * np.exp(-1.6 * r) * (1 + 0.3 * r)
        return u_s, w_d
    V = get_deuteron_potential(waveform, r_vals_fm)
    # V_central = np.interp(r, r_vals_fm, V).reshape((N,N,N))


    # Initial wavefunctions
    psi_S = psi_S_initial.astype(np.complex128)
    psi_D = psi_D_initial.astype(np.complex128)

    # Absorbing boundary mask
    absorbing_mask = get_absorbing_boundary_mask((N,N,N), 10, 0.5)


    # Create the figure with slider
    fig = go.Figure()

    all_S_densities = []
    all_D_densities = []

    for i in range(num_steps):
        psi_S, psi_D = split_step_coupled(psi_S, psi_D, V_central, V_tensor, dt, absorbing_mask)
        if i % plot_every == 0:
            all_S_densities.append(np.abs(psi_S)**2)
            all_D_densities.append(np.abs(psi_D)**2)

    for step_idx, (s_density, d_density) in enumerate(zip(all_S_densities, all_D_densities)):
        iso_level_s = np.mean(s_density) + 2*np.std(s_density)
        iso_level_d = np.mean(d_density) + 2*np.std(d_density)

        fig.add_trace(go.Isosurface(
            x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
            value=s_density.flatten(),
            isomin=iso_level_s, isomax=iso_level_s,
            surface_count=1,
            caps=dict(x_show=False, y_show=False),
            colorscale='Reds',
            name=f'S-wave (t={step_idx*plot_every})',
            visible=(step_idx == 0)
        ))
        fig.add_trace(go.Isosurface(
            x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
            value=d_density.flatten(),
            isomin=iso_level_d, isomax=iso_level_d,
            surface_count=1,
            caps=dict(x_show=False, y_show=False),
            colorscale='Blues',
            name=f'D-wave (t={step_idx*plot_every})',
            visible=(step_idx == 0)
        ))

    # Create slider
    steps = []
    for i in range(len(all_S_densities)):
        step = dict(
            method="update",
            args=[{"visible": [False] * len(fig.data)}],
            label=f"{i * plot_every}"
        )
        step["args"][0]["visible"][i*2] = True
        step["args"][0]["visible"][i*2 + 1] = True
        steps.append(step)

    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Time Step: "},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_layout(
        title="Time Evolution of S- and D-Wave Probability Densities",
        scene=dict(
            xaxis_title="X (fm)",
            yaxis_title="Y (fm)",
            zaxis_title="Z (fm)"
        ),
        sliders=sliders
    )

    fig.show()


