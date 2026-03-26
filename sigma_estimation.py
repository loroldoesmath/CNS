"""
Estimate Wilson-Cowan noise amplitudes (sigma_E, sigma_I) from EEG data.

APPROACH
--------
The linearized stochastic Wilson-Cowan model predicts a temporal power
spectral density (PSD) for the excitatory field X(x,t). At the k=0
(spatially uniform) mode, which dominates scalp EEG, the predicted PSD is:

    S_model(omega) = sigma_E^2 * G_EE(omega) + sigma_I^2 * G_II(omega)

where G_EE and G_II are transfer functions determined entirely by the
Jacobian A(0) of the E-I system (i.e., the E-I parameters alone, no sigma).

Since S_model is *linear* in (sigma_E^2, sigma_I^2), fitting sigma reduces
to a non-negative least squares problem:

    minimize  ||S_observed(omega) - [G_EE(omega), G_II(omega)] [s_E, s_I]^T ||^2
    subject to  s_E, s_I >= 0

where s_E = sigma_E^2, s_I = sigma_I^2.

USAGE
-----
    python estimate_sigma.py path/to/file.edf

DEPENDENCIES
------------
    pip install mne numpy scipy matplotlib
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from scipy.signal import welch


# -------------------------------------------------------------------------
# 1. Wilson-Cowan model parameters
#    These are fixed. Adjust to match your chosen parameter regime.
#    The key requirement: A(0) should have eigenvalues with small negative
#    real part and nonzero imaginary part (stable, near-Hopf regime).
# -------------------------------------------------------------------------

p = {
    'tau_E':  1.0,    # E time constant (arbitrary units)
    'tau_I':  2.0,    # I time constant
    'w_EE':  14.0,    # E->E weight  (near-Hopf: increase to approach bifurcation)
    'w_EI':  12.0,    # I->E weight
    'w_IE':  15.0,    # E->I weight
    'w_II':   3.0,    # I->I weight
    'beta_E':  1.0,   # E sigmoid steepness
    'theta_E': 4.0,   # E sigmoid threshold
    'beta_I':  1.0,   # I sigmoid steepness
    'theta_I': 3.7,   # I sigmoid threshold
    'P': 1.25,        # external input to E
    'Q': 0.0,         # external input to I
}


# -------------------------------------------------------------------------
# 2. Model utilities
# -------------------------------------------------------------------------

def sigmoid(u, beta, theta):
    return 1.0 / (1.0 + np.exp(-beta * (u - theta)))

def sigmoid_prime(u, beta, theta):
    s = sigmoid(u, beta, theta)
    return beta * s * (1.0 - s)

def find_steady_state(p):
    """Find the homogeneous steady state (E*, I*) numerically."""
    from scipy.optimize import fsolve
    def rhs(state):
        E, I = state
        dE = -E + sigmoid(p['w_EE']*E - p['w_EI']*I + p['P'], p['beta_E'], p['theta_E'])
        dI = -I + sigmoid(p['w_IE']*E - p['w_II']*I + p['Q'], p['beta_I'], p['theta_I'])
        return [dE, dI]
    best = None
    best_res = np.inf
    for E0 in np.linspace(0.05, 0.95, 8):
        for I0 in np.linspace(0.05, 0.95, 8):
            sol = fsolve(rhs, [E0, I0], full_output=True)
            res = np.max(np.abs(rhs(sol[0])))
            if sol[2] == 1 and res < best_res:
                best_res = res
                best = sol[0]
    return best

def jacobian_k0(p, E_ss, I_ss):
    """Jacobian of linearized WC system at (E*, I*) for k=0 (no diffusion)."""
    u_E = p['w_EE']*E_ss - p['w_EI']*I_ss + p['P']
    u_I = p['w_IE']*E_ss - p['w_II']*I_ss + p['Q']
    SE = sigmoid_prime(u_E, p['beta_E'], p['theta_E'])
    SI = sigmoid_prime(u_I, p['beta_I'], p['theta_I'])
    A = np.array([
        [(-1 + p['w_EE']*SE) / p['tau_E'],  (-p['w_EI']*SE) / p['tau_E']],
        [( p['w_IE']*SI) / p['tau_I'],       (-1 - p['w_II']*SI) / p['tau_I']]
    ])
    return A


# -------------------------------------------------------------------------
# 3. Analytical PSD from linearized model
#
#    For the k=0 mode, the stochastic WC equations in Fourier time give:
#
#      d/dt [X_hat, v_hat] = A(0) [X_hat, v_hat] + [sigma_E/tau_E, sigma_I/tau_I] noise
#
#    The PSD of X_hat at temporal frequency omega is the (1,1) entry of:
#      S(omega) = (-iw*I - A)^{-1} Q (iw*I - A^T)^{-1}
#
#    Since Q = diag(sigma_E^2/tau_E^2, sigma_I^2/tau_I^2), and S is linear in
#    the diagonal entries of Q, we can factor:
#      S_X(omega) = sigma_E^2 * G_EE(omega) + sigma_I^2 * G_II(omega)
#
#    where G_EE and G_II are computed below with sigma=1.
# -------------------------------------------------------------------------

def transfer_matrix(A, omega):
    """Resolvent (-i*omega*I - A)^{-1} as a 2x2 complex matrix."""
    M = -1j * omega * np.eye(2) - A
    return np.linalg.inv(M)

def psd_components(A, omega_arr, tau_E, tau_I):
    """
    Compute G_EE(omega) and G_II(omega):
    the contribution to S_X(omega) per unit sigma_E^2 and sigma_I^2.

    Returns
    -------
    G_EE : (N_omega,) array
    G_II : (N_omega,) array
    """
    G_EE = np.zeros(len(omega_arr))
    G_II = np.zeros(len(omega_arr))

    # Unit noise covariance matrices (one for each source)
    Q_E = np.array([[1.0/tau_E**2, 0.0], [0.0, 0.0]])
    Q_I = np.array([[0.0, 0.0], [0.0, 1.0/tau_I**2]])

    for i, omega in enumerate(omega_arr):
        R = transfer_matrix(A, omega)          # (-iw I - A)^{-1}
        Rh = transfer_matrix(A, omega).conj().T  # (iw I - A^T)^{-1} = R^H

        S_E = R @ Q_E @ Rh
        S_I = R @ Q_I @ Rh

        # (1,1) entry = contribution to excitatory field PSD
        G_EE[i] = np.real(S_E[0, 0])
        G_II[i] = np.real(S_I[0, 0])

    return G_EE, G_II


# -------------------------------------------------------------------------
# 4. Load EEG from EDF and compute empirical PSD
# -------------------------------------------------------------------------

def load_edf_and_compute_psd(edf_path, duration_s=60.0, nperseg_s=4.0,
                              freq_max_hz=40.0):
    """
    Load an EDF file, pick a representative EEG channel, and compute
    the Welch PSD.

    Parameters
    ----------
    edf_path    : path to .edf file
    duration_s  : how many seconds to use (from the start; set None for all)
    nperseg_s   : Welch segment length in seconds
    freq_max_hz : maximum frequency to retain

    Returns
    -------
    freqs_hz    : (N,) frequency array in Hz
    psd         : (N,) power spectral density (uV^2/Hz), averaged over channels
    sfreq       : sampling frequency
    ch_names    : list of channel names used
    """
    try:
        import mne
    except ImportError:
        raise ImportError("Install mne: pip install mne")

    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)
    sfreq = raw.info['sfreq']

    # Pick EEG channels only
    raw.pick_types(eeg=True, verbose=False)
    if len(raw.ch_names) == 0:
        raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)
        print("No EEG channels found; using all channels.")

    ch_names = raw.ch_names
    print(f"Channels: {ch_names}")
    print(f"Sampling frequency: {sfreq} Hz")
    print(f"Duration: {raw.times[-1]:.1f} s")

    # Crop to requested duration
    if duration_s is not None:
        t_end = min(duration_s, raw.times[-1])
        raw.crop(tmin=0, tmax=t_end)

    data, _ = raw[:]   # (n_channels, n_times) in volts

    # Compute Welch PSD for each channel, then average
    nperseg = int(nperseg_s * sfreq)
    psds = []
    for ch_data in data:
        f, pxx = welch(ch_data, fs=sfreq, nperseg=nperseg)
        psds.append(pxx)
    psd_mean = np.mean(psds, axis=0)  # average over channels

    # Convert V^2/Hz -> uV^2/Hz
    psd_mean *= 1e12

    # Restrict to [0, freq_max_hz]
    mask = f <= freq_max_hz
    return f[mask], psd_mean[mask], sfreq, ch_names


# -------------------------------------------------------------------------
# 5. Match model PSD to observed PSD: estimate sigma_E, sigma_I
#
#    The model PSD is in model time units (tau_E = 1 arbitrary unit).
#    The observed PSD is in Hz (1/seconds).
#    We need a time-unit conversion factor: tau_scale (seconds per model unit).
#    This is an additional free parameter; we grid-search over it.
#
#    For each tau_scale candidate:
#      1. Convert observed frequencies from Hz to model units:
#           omega_model = 2 * pi * f_hz * tau_scale
#      2. Compute G_EE(omega_model) and G_II(omega_model)
#      3. Solve NNLS: [sigma_E^2, sigma_I^2] = nnls([G_EE, G_II], S_obs)
#      4. Record residual
#    Return the best-fit tau_scale and sigma values.
# -------------------------------------------------------------------------

def estimate_sigma(freqs_hz, psd_obs, A, tau_E, tau_I,
                   tau_scale_range=(0.005, 0.5), n_tau=80,
                   freq_min_hz=1.0, freq_max_hz=30.0):
    """
    Estimate sigma_E, sigma_I by fitting the model PSD to observed EEG PSD.

    Parameters
    ----------
    freqs_hz     : observed frequency array (Hz)
    psd_obs      : observed PSD (uV^2/Hz)
    A            : 2x2 Jacobian matrix at k=0
    tau_E, tau_I : model time constants
    tau_scale_range : (min, max) range for tau conversion factor (seconds)
    n_tau        : number of tau_scale candidates to try
    freq_min_hz  : lower frequency bound for fitting
    freq_max_hz  : upper frequency bound for fitting

    Returns
    -------
    sigma_E, sigma_I : estimated noise amplitudes
    tau_scale        : best-fit time-unit conversion (seconds per model unit)
    fit_info         : dict with residuals and arrays for plotting
    """
    # Restrict fitting to [freq_min_hz, freq_max_hz]
    mask = (freqs_hz >= freq_min_hz) & (freqs_hz <= freq_max_hz)
    f_fit = freqs_hz[mask]
    S_fit = psd_obs[mask]

    # Normalize observed PSD for numerical stability
    S_scale = np.max(S_fit)
    S_fit_norm = S_fit / S_scale

    tau_scales = np.linspace(tau_scale_range[0], tau_scale_range[1], n_tau)
    best_res = np.inf
    best = None

    for tau_scale in tau_scales:
        # Convert Hz to model angular frequency units
        omega_model = 2 * np.pi * f_fit * tau_scale

        G_EE, G_II = psd_components(A, omega_model, tau_E, tau_I)

        # Clip negative values (numerical artifact far from resonance)
        G_EE = np.clip(G_EE, 0, None)
        G_II = np.clip(G_II, 0, None)

        # Build design matrix: columns are G_EE and G_II
        B = np.column_stack([G_EE, G_II])

        # Normalize columns
        col_scales = np.maximum(np.max(np.abs(B), axis=0), 1e-30)
        B_norm = B / col_scales

        # Non-negative least squares: fit [s_E, s_I] = [sigma_E^2, sigma_I^2]
        coeffs_norm, residual = nnls(B_norm, S_fit_norm)
        coeffs = coeffs_norm / col_scales  # undo column scaling

        # Residual as fraction of signal power
        S_pred = B @ coeffs
        rel_res = np.sum((S_fit_norm - B_norm @ coeffs_norm)**2) / np.sum(S_fit_norm**2)

        if rel_res < best_res:
            best_res = rel_res
            s_E, s_I = coeffs * S_scale  # restore physical scale
            best = {
                'sigma_E': np.sqrt(max(s_E, 0)),
                'sigma_I': np.sqrt(max(s_I, 0)),
                'tau_scale': tau_scale,
                'rel_residual': rel_res,
                'G_EE': G_EE,
                'G_II': G_II,
                'S_pred': S_pred * S_scale,
                'f_fit': f_fit,
                'S_fit': S_fit,
            }

    return best


# -------------------------------------------------------------------------
# 6. Plotting
# -------------------------------------------------------------------------

def plot_results(freqs_hz, psd_obs, fit_info, A, eig_vals):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle('Wilson-Cowan Noise Amplitude Estimation from EEG', fontsize=12)

    # Panel 1: observed vs fitted PSD
    ax = axes[0]
    ax.semilogy(freqs_hz, psd_obs, color='steelblue', lw=1, alpha=0.7, label='Observed EEG PSD')
    ax.semilogy(fit_info['f_fit'], fit_info['S_pred'], color='tomato', lw=2,
                label=f"Model fit\n$\\sigma_E$={fit_info['sigma_E']:.3f}, "
                      f"$\\sigma_I$={fit_info['sigma_I']:.3f}")
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (μV²/Hz)')
    ax.set_title('Observed vs. fitted PSD')
    ax.legend(fontsize=9)

    # Panel 2: transfer functions G_EE and G_II
    ax = axes[1]
    ax.plot(fit_info['f_fit'], fit_info['G_EE'], label='$G_{EE}$ (E noise → E field)', color='coral')
    ax.plot(fit_info['f_fit'], fit_info['G_II'], label='$G_{II}$ (I noise → E field)', color='mediumpurple')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Transfer function (a.u.)')
    ax.set_title('Model transfer functions')
    ax.legend(fontsize=9)

    # Panel 3: eigenvalues of A(0) in complex plane
    ax = axes[2]
    ax.axvline(0, color='k', lw=0.8, ls='--', alpha=0.5)
    ax.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
    for ev in eig_vals:
        ax.plot(ev.real, ev.imag, 'o', ms=10, color='steelblue')
        ax.annotate(f'  λ={ev:.3f}', (ev.real, ev.imag), fontsize=9)
    ax.set_xlabel('Re(λ)')
    ax.set_ylabel('Im(λ)')
    ax.set_title('Eigenvalues of A(0)\n(stability: Re < 0 required)')
    ax.set_xlim([min(ev.real for ev in eig_vals) * 1.5, 0.3])

    plt.tight_layout()
    plt.savefig('sigma_estimation.png', dpi=150)
    plt.show()
    print("Saved sigma_estimation.png")


# -------------------------------------------------------------------------
# 7. Main
# -------------------------------------------------------------------------

def main(edf_path):
    print("=" * 60)
    print("Wilson-Cowan sigma estimation from EEG")
    print("=" * 60)

    # Step 1: Find steady state and Jacobian
    print("\n[1] Computing steady state and Jacobian...")
    E_ss, I_ss = find_steady_state(p)
    print(f"    Steady state: E* = {E_ss:.4f}, I* = {I_ss:.4f}")

    A = jacobian_k0(p, E_ss, I_ss)
    eig_vals = np.linalg.eigvals(A)
    print(f"    Eigenvalues of A(0):")
    for ev in eig_vals:
        stability = "STABLE" if ev.real < 0 else "UNSTABLE"
        print(f"      λ = {ev:.4f}  [{stability}]")

    if any(ev.real >= 0 for ev in eig_vals):
        print("\n    WARNING: fixed point is unstable. Reduce w_EE to enter stable regime.")
        print("    Linearized sigma estimation is only valid in the stable regime.")

    # Step 2: Load EEG and compute PSD
    print(f"\n[2] Loading EEG from {edf_path}...")
    freqs_hz, psd_obs, sfreq, ch_names = load_edf_and_compute_psd(
        edf_path, duration_s=60.0, nperseg_s=4.0, freq_max_hz=40.0
    )
    print(f"    PSD computed over {len(freqs_hz)} frequency bins up to {freqs_hz.max():.1f} Hz")

    # Step 3: Estimate sigma
    print("\n[3] Fitting model PSD to observed PSD...")
    print("    (grid search over time-unit conversion factor tau_scale)")
    fit = estimate_sigma(
        freqs_hz, psd_obs, A, p['tau_E'], p['tau_I'],
        tau_scale_range=(0.005, 0.5), n_tau=100,
        freq_min_hz=1.0, freq_max_hz=30.0
    )

    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"  sigma_E       = {fit['sigma_E']:.4f}  (excitatory noise amplitude)")
    print(f"  sigma_I       = {fit['sigma_I']:.4f}  (inhibitory noise amplitude)")
    print(f"  tau_scale     = {fit['tau_scale']:.4f} s  (model time unit in seconds)")
    print(f"  Relative residual = {fit['rel_residual']:.4f}")
    print()
    print("  Interpretation:")
    print(f"    One model time unit ≈ {fit['tau_scale']*1000:.1f} ms")
    dominant = 'excitatory' if fit['sigma_E'] > fit['sigma_I'] else 'inhibitory'
    print(f"    Dominant noise source: {dominant} population")
    print()
    print("  Caveats:")
    print("    - E-I parameters are fixed at assumed values; sigma estimates")
    print("      are conditional on this choice.")
    print("    - tau_scale is a free parameter absorbing unit ambiguity.")
    print("    - The fit uses the k=0 (spatially uniform) mode only, consistent")
    print("      with EEG as a spatial average.")

    # Step 4: Plot
    print("\n[4] Plotting...")
    plot_results(freqs_hz, psd_obs, fit, A, eig_vals)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python estimate_sigma.py path/to/file.edf")
        print("\nRunning on synthetic EEG data for demonstration...")

        # --- Synthetic demo: generate fake EEG-like PSD and run estimation ---
        np.random.seed(42)
        freqs_hz = np.linspace(0.5, 40, 300)
        # Fake PSD: 1/f background + alpha peak at ~10 Hz
        psd_fake = (100.0 / (freqs_hz + 1)**1.5
                    + 20.0 * np.exp(-0.5 * ((freqs_hz - 10) / 2)**2)
                    + 2.0 * np.random.exponential(1, len(freqs_hz)))

        E_ss, I_ss = find_steady_state(p)
        A = jacobian_k0(p, E_ss, I_ss)
        eig_vals = np.linalg.eigvals(A)

        print(f"\nSteady state: E* = {E_ss:.4f}, I* = {I_ss:.4f}")
        print(f"Eigenvalues: {eig_vals}")

        fit = estimate_sigma(freqs_hz, psd_fake, A, p['tau_E'], p['tau_I'])

        print(f"\nsigma_E = {fit['sigma_E']:.4f}")
        print(f"sigma_I = {fit['sigma_I']:.4f}")
        print(f"tau_scale = {fit['tau_scale']:.4f} s")
        print(f"Relative residual = {fit['rel_residual']:.4f}")

        plot_results(freqs_hz, psd_fake, fit, A, eig_vals)
    else:
        main(sys.argv[1])