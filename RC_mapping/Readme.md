# RC Mapping: Continuously Monitored Non-Markovian Quantum System

This folder contains a Jupyter notebook implementation of the RC (reaction coordinate) mapping analysis used in the manuscript:

**“Parameter Estimation in a Continuously Monitored Non-Markovian Quantum System.”**

The main notebook `RC_mapping.ipynb` compares three dynamical descriptions for a harmonic system coupled to a structured environment:

1. **LOCAL**: GKLS dynamics for the system mode only  
2. **GLOBAL**: GKLS dynamics for the enlarged system + reaction coordinate (S-RC unit)  
3. **EXACT**: Non-Markovian dynamics from inverse Laplace methods and noise-kernel integration

---

## Model and Parameters

The notebook initializes (example values shown in the code):

- Bare system frequency: `Omega_S_bare = 1.0`
- Environment parameters: `alpha1`, `alpha2`, `gamma`
- Coupling and RC frequency:
  - `lam = sqrt(alpha1 * alpha2 * gamma)`
  - `Omega_RC = sqrt(alpha2 * gamma)`
- Effective system renormalization:
  - `dOmega_S = lam / Omega_RC`
  - `Omega_S = sqrt(Omega_S_bare^2 + dOmega_S^2)`
- Temperature: `temp`
- Time grid: `t in [0, 50]` with `dt = 1e-3`

---

## Section Overview

## 1) LOCAL (GKLS for S)

Implements a 2D Gaussian dynamics (`x, p`) with drift and diffusion matrices:

- `lindblad_matrices_system(...)`
- `lindblad_ode_system(...)`
- `solve_system(...)`

The damping rate is derived from the structured spectral density:
- `spectral_density_original(omega, gamma, lam, Omega_RC)`
- `kappa = J(Omega_S_bare)/Omega_S_bare`

Outputs:
- `mean_t_local`
- `cov_t_local`
- `var_xS_local`, `var_pS_local`, `xS_pS_local`

Also computes steady-state covariance with:
- `solve_continuous_lyapunov(A, -D)`

---

## 2) GLOBAL (GKLS for S-RC)

Implements 4D Gaussian dynamics (`x_S, p_S, x_RC, p_RC`) for the enlarged Markovian embedding:

- normal-mode frequencies via `frequencies(...)`
- drift/diffusion via `lindblad_matrices_system_RC(...)`
- time evolution via `solve_system_RC(...)` using `solve_ivp(..., method="Radau")`

Outputs:
- `mean_t`, `cov_t`
- `var_xS_global`, `var_pS_global`
- `var_xRC_global`, `var_pRC_global`
- `xSpS_global`

Steady-state covariance is again extracted from:
- `solve_continuous_lyapunov(A, -D)`

---

## 3) EXACT (non-Markovian reference)

Constructs the Green function through partial fraction expansion of a rational Laplace transform:

- `inverse_laplace_rational(B, A)`
- `derivatives_from_residues(...)` -> `G(t), Gdot(t), Gddot(t)`

Computes the thermal noise kernel using residue + Matsubara contributions:

- `noise_kernel_residue(t, gamma, alpha1, alpha2, Temp, Nmats=...)`

Then propagates means and covariances by combining:
- deterministic Green-function part
- bath-induced double-convolution terms (Toeplitz acceleration via `matmul_toeplitz`)

Main routine:
- `evolve_dynamics(...)`

Outputs:
- `sol["x_mean"]`, `sol["p_mean"]`
- `sol["var_x"]`, `sol["var_p"]`
- `sol["x_p"]`

---

## Plots Produced

The notebook generates comparison plots for:

- `⟨x⟩`: global GKLS vs local GKLS vs exact
- `⟨p⟩`: global GKLS vs local GKLS vs exact
- `⟨x²⟩`: global GKLS vs local GKLS vs exact
- `⟨p²⟩`: global GKLS vs local GKLS vs exact
- `⟨xp⟩`: global GKLS vs local GKLS vs exact

Styling uses `scienceplots` and LaTeX text rendering.

---

## Requirements

Typical Python dependencies used in the notebook:

- `numpy`
- `scipy`
- `matplotlib`
- `mpmath`
- `tqdm`
- `scienceplots` (for final figure style)
