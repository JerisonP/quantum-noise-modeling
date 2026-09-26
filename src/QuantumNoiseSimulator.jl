"""
    QuantumNoiseSimulator

Generate and validate classical noise processes ζ(t) that drive a qubit
Hamiltonian, and (in later modules) propagate the qubit under that noise.

# Layout

| File                        | Contents                                               |
|-----------------------------|--------------------------------------------------------|
| `noise/ensemble.jl`         | `NoiseEnsemble` — the container every generator returns |
| `noise/interface.jl`        | `AbstractNoiseModel` and the functions every model implements |
| `noise/ou.jl`               | `OUNoiseModel` — Ornstein–Uhlenbeck, exact discretisation |
| `noise/white.jl`            | `WhiteNoiseModel` — i.i.d. Gaussian samples             |
| `noise/fractional.jl`       | `FractionalNoiseModel` — 1/f^α noise (Kasdin 1995)      |

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

using FFTW: plan_rfft, plan_irfft
using Random: AbstractRNG, default_rng
using SpecialFunctions: gamma

include("noise/ensemble.jl")
include("noise/interface.jl")
include("noise/ou.jl")
include("noise/white.jl")
include("noise/fractional.jl")

# Container
export NoiseEnsemble, times, samples, ntimes, ntrajectories, timestep, subensemble

# Model interface
export AbstractNoiseModel, generate_ensemble
export noise_mean, stationary_variance, theoretical_autocovariance, theoretical_psd

# Models
export OUNoiseModel, WhiteNoiseModel, FractionalNoiseModel

# 1/f^α building blocks (used by the validation layer and the tests)
export driving_variance, pulse_response, ar_coefficients, fir_filter, ar_filter

end # module
