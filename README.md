# QuantumNoiseSimulator

[![CI](https://github.com/OWNER/REPO/actions/workflows/CI.yml/badge.svg)](https://github.com/OWNER/REPO/actions/workflows/CI.yml)

Generate classical noise ζ(t) for a driven qubit, **validate that it has the
statistics it claims**, and propagate the qubit under it.

Noise models: Ornstein–Uhlenbeck (exact discretisation), white, 1/f^α (Kasdin
fractional differencing, including the nonstationary pink-noise case α = 1), and
band-limited stationary 1/f (matching the Mathematica kernel Ci(100|τ|) − Ci(|τ|)).

> 🚧 This is a ground-up rewrite, in progress. See [`docs/ROADMAP.md`](docs/ROADMAP.md)
> for the step plan.

## Install

```julia
julia> ]
pkg> activate .
pkg> instantiate
```

## Quick start

```julia
using QuantumNoiseSimulator, Random

rng  = Xoshiro(2026)                         # explicit RNG ⇒ reproducible
pink = FractionalNoiseModel(1.0, 1e-3)       # α = 1, S(f) ≈ 2·Q_psd/(2πf)
ens  = generate_ensemble(pink, (0.0, 4.0), 1e-3, 1000; rng)

X = samples(ens)          # 4001 × 1000: column j = trajectory j
X[end, :]                 # 1000 i.i.d. draws of ζ(T)

ou = OUNoiseModel(1.0, 0.5)                  # σ = 1, τ_c = 0.5
theoretical_autocovariance(ou, 0.5)          # σ² e^{-1}
theoretical_psd(ou, 1.0)                     # Lorentzian (continuous time)
theoretical_psd(ou, 1.0; dt = 0.01)          # exact PSD of the sampled sequence
```

## Conventions

| Quantity | Convention |
|----------|------------|
| C(τ) | autocovariance E[(ζ(t)−μ)(ζ(t+τ)−μ)]; ρ = C/C(0) is the autocorrelation |
| S(f) | **one-sided**, f in cycles per unit time, ∫₀^∞ S df = variance |
| 1/f^α amplitude | `Q_psd` is Kasdin's Q: S(f) ≈ 2·Q_psd/(2πf)^α (one-sided); innovations Q_d = Q_psd·Δt^(α−1) |
| `dt` in theory functions | `nothing` → continuous process; a value → exactly the sampled sequence |
| Randomness | every generator takes `rng`; tests use `StableRNGs` |

## Repository layout

```
src/
  QuantumNoiseSimulator.jl   module, exports
  noise/                     ensemble container, model interface, OU / white / 1/f^α / band-limited 1/f
  analysis/                  estimators, fits, exact expectations, validate_noise_model
test/                        one test file per source file + Aqua hygiene checks
docs/                        explainer per step, roadmap, REFERENCES (formula → paper)
notebooks/                   validation appendix (rebuilt in Step 6)
```

## Running the tests

```julia
pkg> test
```

CI runs the same suite on Julia 1.10 (LTS) and the latest release on every
push and pull request.

## References

- N. J. Kasdin, "Discrete simulation of colored noise and stochastic processes
  and 1/f^α power law noise generation", *Proc. IEEE* **83**, 802 (1995).
- J. R. M. Hosking, "Fractional differencing", *Biometrika* **68**, 165 (1981).
- P. Dutta and P. M. Horn, "Low-frequency fluctuations in solids: 1/f noise",
  *Rev. Mod. Phys.* **53**, 497 (1981).
- J. Ruseckas and B. Kaulakys, "1/f noise from nonlinear stochastic differential
  equations", *Phys. Rev. E* **81**, 031105 (2010).
- Wikipedia, "Ornstein–Uhlenbeck process".

See [`docs/REFERENCES.md`](docs/REFERENCES.md) for which equation each line of code implements.