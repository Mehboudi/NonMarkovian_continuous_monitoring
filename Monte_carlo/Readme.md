# Monte Carlo folder

Run the code from:

- `main_runner.m`

This is the entry point. Most files in this folder are called by it (directly or indirectly).

## What this code does

`main_runner.m` runs Bayesian temperature estimation under continuous monitoring using Monte Carlo trajectories.  
The script:

1. sets simulation/prior controls,
2. initializes physical and derived quantities,
3. precomputes BCRB/Fisher-information-related quantities,
4. runs one tracked trajectory + many Monte Carlo runs,
5. plots and saves outputs.

Core calls from `main_runner.m`:
- `init_physical_parameters`
- `compute_derived_quantities`
- `compute_constant_matrices`
- `compute_bcrb`
- `simulate_true_system`
- `update_bayesian_filter`
- `plot_simulation_results`

---

## Where to set parameters

## 1) Main run controls (first place to edit)
In `main_runner.m` (top block “Simulation Controls”), set:

- `estimation_metric` (`'EMSLE'` or `'ERMSE'`)
- `T_true_nominal`
- `t_final`
- `dt`
- `N_hypotheses`
- `num_m`
- prior range/shape: `T_min`, `T_max`, `prior_alpha`

These strongly affect runtime and memory.  
Large `t_final/dt` and large `num_m` produce very large outputs.

## 2) Physical constants
In `init_physical_parameters.m`:
- `eps`
- `alpha1`
- `alpha2`
- `gamma`
- `Lambda_cutoff`
- `lambda_meas` (measurement strength)

## 3) Derived model quantities
In `compute_derived_quantities.m`:
- base frequency: `w_s_bare`
- derived frequencies/mixing are built here (`w_rc`, `w_p`, `w_m`, `theta`, etc.)
- bosonic occupation and derivative functions are defined here:
  - `N_bosonic(w,T)`
  - `dN_dT(w,T)`

## 4) Constant/derived matrices (including measurement + environment covariance usage)
In `compute_constant_matrices.m` and temperature-dependent builder `evaluate_T_matrices.m`:

- `evaluate_T_matrices.m` builds temperature-dependent diffusion/process-noise pieces:
  - `D(T)`, `D'(T)`, `Q(T)`.

- In `compute_bcrb.m`, the environment covariance block is explicitly set as
  `Sigma_env_full(3:4,3:4) = 0.5*eye(2)` and used in the steady-state Riccati setup.

- Measurement-related matrices (`C_meas_full`, `V_inv`, `V_inv_sqrt`) are read from `sys_matrices` and used in:
  - `compute_bcrb.m`,
  - `simulate_true_system.m`,
  - `update_bayesian_filter.m`.

So if you want to change measurement model / covariance structure, check matrix construction in `compute_constant_matrices.m` and related use in those files.

## 5) Initial state/covariance for true system and hypotheses
In `main_runner.m`:
- true system initialization:
  - `true_sys.d = zeros(4,1)`
  - `true_sys.sigma = 0.5 * diag([coth(...), coth(...), 1, 1])` pattern
- hypothesis initialization:
  - `hypotheses.d = zeros(4, N_hypotheses)`
  - `hypotheses.sigma(:,:,i1) = ...` (temperature-dependent diagonal initialization)
  - `hypotheses.weights` initialized from prior.

If you need different initial matrices, this is one of the key places to edit.

---

## Output files and size warning

`main_runner.m` saves:
- `matlab.mat` (`-v7.3`)
- `matlab_plot_payload.mat` (`-v7.3`)

The data can become huge (especially with large `N_steps`, `num_m`, or dense saved histories).

---

## Prior / prior-information utilities (including files not directly in main run)

Main pipeline prior utilities:
- `prior_utils.m`  
  (prior model creation, prior PDF, derivative, sampling, normalization)
- used directly by `main_runner.m` and `compute_bcrb.m`.

Standalone prior analysis script not called by `main_runner.m`:
- `Prior_plot.m`  
  (independent script to visualize the prior distribution).

---

## Practical way to adapt to a new problem

Start from `main_runner.m` and follow calls in order.  
For most adaptations:

1. Edit run controls in `main_runner.m`.
2. Edit physical constants in `init_physical_parameters.m`.
3. Edit matrix definitions in `compute_constant_matrices.m` / `evaluate_T_matrices.m`.
4. If needed, edit initialization of `true_sys` and `hypotheses` in `main_runner.m`.
