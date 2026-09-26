# Step 1: Noise models

*What this layer does, the math behind each generator, and what the tests prove.*

## 1. The job of this layer

The qubit feels a random field ζ(t). Every question downstream (spectra, gate
fidelity) depends only on the **statistics** of ζ: its mean, its
autocovariance C(τ), or equivalently its power spectral density S(f). So a
noise model is a recipe for producing sample paths with *known* statistics.
"Validating" a generator means showing that the samples really have those statistics.

Every model answers two questions:

| Question | Function |
|----------|----------|
| "Give me M sample paths on this time grid" | `generate_ensemble(model, (t0, T), dt, M; rng)` |
| "What should their statistics be?" | `noise_mean`, `stationary_variance`, `theoretical_autocovariance`, `theoretical_psd` |

### Conventions (used everywhere)

* **Autocovariance** C(τ) = E[(ζ(t)−μ)(ζ(t+τ)−μ)] has units of ζ². The
  *autocorrelation* is the normalised ρ(τ) = C(τ)/C(0). (The old code called
  C "autocorrelation", which is a common but confusing usage.)
* **One-sided PSD**: f ≥ 0 and ∫₀^∞ S(f) df = variance. The one-sided
  PSD is twice the two-sided one.
* **`dt = nothing` vs `dt` given.** Without `dt`, theory describes the continuous
  process ζ(t). With `dt`, it describes *exactly the sampled sequence* ζ_k the
  generator returns, aliasing included, and the PSD then integrates to the
  variance over [0, 1/(2dt)].
* **Randomness is explicit.** Pass `rng = Xoshiro(seed)` (or a `StableRNG` in
  tests). No global `Random.seed!`, and no dependence on the thread count.

### The ensemble container

`NoiseEnsemble(t, X)` holds a uniform time grid `t` (length n_t) and an
`n_t × M` matrix `X` whose **column j is trajectory j**. Two slices do most of the work:

* `X[:, j]`: one path in time, used for time averages, periodograms, and ACFs.
* `X[k, :]`: M *independent* draws of ζ(t_k), used for ensemble averages at a fixed time.

The second slice matters because its entries are i.i.d., so textbook
standard errors (σ/√M for a mean, σ²√(2/(M−1)) for a variance) are *exact*.
The tests lean on this.

## 2. Ornstein–Uhlenbeck (OU)

    dζ = −(ζ − μ)/τ_c dt + √(2σ²/τ_c) dW

One correlation time τ_c. The stationary law is N(μ, σ²),
C(τ) = σ² e^{−|τ|/τ_c}, and the PSD is a Lorentzian 4σ²τ_c/(1+(2πfτ_c)²): flat below
f_c = 1/(2πτ_c) and falling as f⁻² above it.

**How we generate it, exactly.** The OU process is Gaussian and Markov, and its
transition over a step Δt is known in closed form:

    ζ_{k+1} − μ = ϕ(ζ_k − μ) + σ√(1−ϕ²) ξ_k,     ϕ = e^{−Δt/τ_c},  ξ_k ~ N(0,1).

Start from ζ_0 ~ N(μ, σ²), which is already stationary, and iterate. This is
**exact**: the samples have precisely the OU joint distribution for any Δt,
with no solver and no tolerance. The old code solved the SDE with an adaptive integrator,
which is slower, not reproducible across thread counts, and has a small solver error.
Because the update is exact, the old "time-step convergence" section is replaced by
tests that the recursion reproduces ϕ exactly.

**Sampled spectrum.** The PSD of the sampled sequence is not exactly the Lorentzian: power
above Nyquist folds back (aliasing). Its exact form is
2Δtσ²(1−ϕ²)/(1−2ϕcos 2πfΔt + ϕ²), available as `theoretical_psd(m, f; dt)`. It
integrates to σ² over [0, 1/2Δt] and approaches the Lorentzian for f ≪ 1/2Δt,
and both properties are tested.

## 3. White noise

ζ_k ~ N(μ, σ²), i.i.d. The PSD of the sequence is flat at 2σ²Δt. The strength
therefore depends on Δt: halving Δt at fixed σ halves the PSD. To keep the physical
PSD fixed, scale σ ∝ 1/√Δt. Continuous white noise has infinite variance, so
`theoretical_psd` requires `dt`.

## 4. 1/f^α noise (Kasdin fractional differencing)

### The construction
Feed i.i.d. innovations w_k ~ N(0, Q_d) through the filter
H(z) = (1 − z⁻¹)^{−α/2}:

* α = 0 gives H = 1, so the output is white.
* α = 2 gives H = 1/(1 − z⁻¹), a running sum, so the output is a random walk (Brownian, 1/f²).
* α = 1 is the geometric midpoint: pink, 1/f.

The output PSD is |H|² × (innovation PSD):

    S(f) = 2 Q_d Δt / |2 sin(πfΔt)|^α   →   2 Q_psd / (2πf)^α   for f ≪ 1/(2Δt).

Choosing **Q_d = Q_psd·Δt^{α−1}** makes the low-frequency level independent of Δt,
so `Q_psd` is the physical amplitude. At α = 1, Q_d = Q_psd whatever Δt is.

### Two ways to apply the filter (both implemented, tested equal)
* **FIR**: x_n = Σ_{k=0}^{n} h_k w_{n−k}, with h_0 = 1 and h_k = h_{k−1}(α/2+k−1)/k.
  This is a convolution, done by zero-padded FFT in O(N log N). The zero padding matters:
  without it the FFT computes a *circular* convolution and the end of the record
  wraps onto the start.
* **AR**: x_n = w_n − Σ_{k≥1} a_k x_{n−k}, with a_0 = 1 and a_k = a_{k−1}(k−1−α/2)/k.
  This is O(N²) and is kept as an independent cross-check.

h and a are inverse power series (H·H⁻¹ = 1), so on the same innovations the two
must agree to rounding error. That is a strong test: an indexing slip in either one
breaks it.

### Stationary or not?
Both filters start from rest (x = 0 for t < 0), so

    Var[x_n] = Q_d Σ_{k=0}^{n} h_k².

* **α < 1**: h_k² ~ k^{α−2} is summable, so the variance converges to
  Q_d Γ(1−α)/Γ(1−α/2)² and the process becomes stationary. `stationary_variance`
  and `theoretical_autocovariance` return these limits. The ACF uses the recursion
  R(k) = R(k−1)(k−1+α/2)/(k−α/2) instead of Kasdin's Γ-function formula, which
  overflows and returns NaN from lag ≈ 172.
* **α ≥ 1**: the sum diverges (∝ ln n at α = 1), so the process is **nonstationary**.
  No stationary C(τ) exists, and the theory functions *throw* instead of returning
  something misleading. What *is* exact for every α is the finite-record prediction,
  which depends on how you estimate. That is Step 3, and it is the fix for
  the notebook's 0.24 ACF RMSE.

## 5. What the tests establish (`test/test_*.jl`)

| Claim | Test |
|-------|------|
| Container rejects malformed data (wrong shape, non-uniform grid, …) | `test_ensemble.jl` |
| Same seed ⇒ identical samples; different seed ⇒ different | OU, fractional |
| OU samples: exact mean/variance at t = 0, 1.5, 3 (stationary from the start) | 5-SE checks on i.i.d. cross-sections |
| OU samples: two-time covariance at lag τ_c is σ²/e; lag-1 regression equals ϕ | 5-SE checks |
| OU/white sampled PSD integrates to σ² (Parseval) and → Lorentzian at low f | deterministic, rtol 1e-6 |
| Kasdin coefficients: known values; α=0 → identity; α=2 → running sum; h·a = 1 | deterministic |
| FIR = naive convolution = AR recursion on the same innovations | rtol 1e-9–1e-12 |
| Generated variance at t = T equals Q_dΣh² for α ∈ {0, 0.5, 1, 1.5, 2} | 5-SE, M = 4000 |
| Stationary theory: recursion = Eq. 110 (k ≤ 50), stays finite at k = 800, ∫S = σ², Σhh → R(m) | deterministic |
| Discrete PSD → continuous PSD as fΔt → 0 | rtol 1e-6 |
| Package hygiene: no ambiguities, stale deps, missing compat, or piracy | Aqua.jl |

**Why 5 standard errors, and what that cannot catch.** A correct generator
fails a given 5-SE check with probability about 6×10⁻⁷, so the suite is not flaky.
The price is resolution: at M = 4000 the relative SE of a variance is
√(2/M) ≈ 2.2%, so only variance errors above roughly 11% are *guaranteed* to fail.
Small errors are caught by the **deterministic** tests instead (FIR ≡ AR,
Parseval integrals, known coefficients), which are exact to rounding. That split
is deliberate: statistical tests catch gross mistakes such as a missing √2 or the
wrong τ, and exact tests catch subtle ones. Before committing, each statistical
check was run over 200 seeds in the NumPy mirror: the z-scores had mean ≈ 0, SD ≈ 1,
and max |z| ≈ 3.4, so the thresholds are calibrated.
