"""
    QuantumNoiseSimulator

Generate and validate classical noise ξ(t), and solve the driven qubit under it by
brute force (thesis Eq. 3.12): H(t) = f_x(t)/2·σx + ξ(t)/2·σz.

# Layout

| File                        | Contents                                               |
|-----------------------------|--------------------------------------------------------|
| `noise/ensemble.jl`         | `NoiseEnsemble` — the container every generator returns |
| `noise/interface.jl`        | `AbstractNoiseModel` and the functions every model implements |
| `noise/ou.jl`               | `OUNoiseModel` — Ornstein–Uhlenbeck, exact discretisation |
| `noise/white.jl`            | `WhiteNoiseModel` — i.i.d. Gaussian samples             |
| `noise/fractional.jl`       | `FractionalNoiseModel` — 1/f^α noise (Kasdin 1995)      |
| `noise/bandlimited.jl`      | `BandLimitedOneOverFNoiseModel` — stationary 1/f in [fl, fh] |
| `analysis/estimators.jl`    | `autocovariance`, `periodogram`, `wiener_khinchin_psd`  |
| `analysis/fits.jl`          | `ols`, `powerlaw_fit`, `exponential_fit`, `student_t_quantile` |
| `analysis/expectations.jl`  | `expected_autocovariance`, `expected_periodogram`, `expected_variance` |
| `analysis/validation.jl`    | `validate_noise_model`, `ValidationReport`                |
| `quantum/qubit.jl`          | Pauli matrices, the six cardinal states, vec/unvec, observables |
| `quantum/gate.jl`           | `CosineGate` and the model H(t) = f_x/2·σx + ξ/2·σz (Eq. 3.12) |
| `quantum/solvers.jl`        | two independent brute-force solvers, `simulate_gate`, `timestep_convergence` |
| `quantum/results.jl`        | ⟨ρ_j(t_g)⟩, eigenvalues, fidelity (Eqs. 3.21, 3.24, 3.25) with standard errors |
| `quantum/small_parameter.jl`| δ (Eq. 3.20) and `delta_sweep`, the reference curves   |
| `quantum/master_equation.jl`| `tcl2_evolution`, the 2nd-order master equation being validated |
| `io/csv.jl`                 | exact CSV save/load: ensembles, channels, δ sweeps (also the hand-in format) |
| `plotting.jl`               | plot function stubs; the methods live in `ext/QuantumNoiseSimulatorPlotsExt.jl` |

# Quick start

```julia
using QuantumNoiseSimulator, Random

rng   = Xoshiro(42)                          # explicit RNG ⇒ reproducible
model = FractionalNoiseModel(1.0, 1e-3)      # α = 1 (pink), Q_psd = 1e-3
ens   = generate_ensemble(model, (0.0, 4.0), 1e-3, 1000; rng)

samples(ens)                                 # 4001 × 1000 matrix, one column per trajectory
```
"""
module QuantumNoiseSimulator

using FFTW: plan_rfft, plan_irfft, rfft
using LinearAlgebra: Hermitian, eigvals, kron, tr
using Random: AbstractRNG, default_rng
using SpecialFunctions: gamma, cosint, beta_inc_inv, erfcinv
using Statistics: mean, std

include("noise/ensemble.jl")
include("noise/interface.jl")
include("noise/ou.jl")
include("noise/white.jl")
include("noise/fractional.jl")
include("noise/bandlimited.jl")
include("analysis/estimators.jl")
include("analysis/fits.jl")
include("analysis/expectations.jl")
include("analysis/validation.jl")
include("quantum/qubit.jl")
include("quantum/gate.jl")
include("quantum/solvers.jl")
include("quantum/results.jl")
include("quantum/small_parameter.jl")
include("quantum/master_equation.jl")
include("io/csv.jl")
include("plotting.jl")

# Container
export NoiseEnsemble, times, samples, ntimes, ntrajectories, timestep, subensemble

# Model interface
export AbstractNoiseModel, generate_ensemble
export noise_mean, stationary_variance, theoretical_autocovariance, theoretical_psd

# Models
export OUNoiseModel, WhiteNoiseModel, FractionalNoiseModel, BandLimitedOneOverFNoiseModel

# 1/f^α building blocks (used by the validation layer and the tests)
export driving_variance, pulse_response, ar_coefficients, fir_filter, ar_filter

# Estimators and fits
export autocovariance, periodogram, wiener_khinchin_psd
export ols, student_t_quantile, powerlaw_fit, exponential_fit

# Exact expectations and validation
export expected_autocovariance, expected_periodogram, expected_variance
export validate_noise_model, ValidationReport, ValidationCheck, passed

# Qubit
export σx, σy, σz, I2, vec_dm, unvec_dm, superoperator
export cardinal_states, bloch_state
export population_0, population_1, purity, bloch_vector
export density_matrix_eigenvalues, validate_density_matrix

# The model (thesis Eq. 3.12)
export CosineGate, drive, rotation_angle, ideal_propagator, ideal_gate
export hamiltonian, noise_axis, interaction_noise

# Brute-force solvers
export su2_exp, propagate_trajectory, simulate_gate, GateResult, timestep_convergence

# Results with standard errors (Eqs. 3.21, 3.24, 3.25)
export apply_channel, final_channel, evolve_state, final_state, final_state_eigenvalues
export fidelity_map, fidelity_trace_formula, trajectory_errors, average_error, average_fidelity
export batch_estimate, choi_matrix, is_cptp

# Small parameter δ (Eq. 3.20) and the 2nd-order master equation under test
export small_noise_parameter, sigma_for_delta, delta_sweep, tcl2_evolution, tcl2_delta_sweep

# Input/output (CSV)
export read_table, save_ensemble, load_ensemble, save_channel, load_channel
export save_delta_sweep, load_delta_sweep

# Plotting (methods are added by the Plots extension)
export plot_noise_traces, plot_validation, plot_state_evolution, plot_timestep_convergence
export plot_error_vs_delta, plot_eigenvalues_vs_delta

end # module