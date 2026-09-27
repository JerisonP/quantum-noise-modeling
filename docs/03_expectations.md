# Step 3: Exact expectations and validation

*The "answer key" for every estimator, and the one function that uses it.*

## 1. The problem this solves

Step 2's estimators are recipes applied to a **finite** record. Their average over
many ensembles (their *expectation*) is not the idealised quantity they are named
after. The old notebook compared estimates with idealised formulas and got
three results that contradicted its own text:

| Old result | Measured vs idealised theory | Measured vs exact expectation |
|------------|-----------------------------|-------------------------------|
| §4.3 pink-noise ACF (normalised RMSE) | 0.24 (text claimed < 0.05) | **0.004** |
| §5.6 OU C(τ_c) with per-record mean removed | 0.27 vs 0.37 | agrees within 1–3 SE |
| §4.4 fitted slope β̂ at α = 0.8 and α = 1.5 | 95% CIs miss α | targets are **0.8025 and 1.5332**, inside the CIs |

The generators were right all along. The comparison target was wrong.

## 2. Everything comes from one matrix

For a record x₀ … x_{N−1}, let Σ[n,m] = Cov(xₙ, xₘ). Every estimator here is
linear or quadratic in x, so its expectation is a simple function of Σ.

| Model | Σ | Why |
|-------|---|-----|
| OU, white, band-limited | **Toeplitz**: Σ[n,m] = c(\|n−m\|) | stationary from t = 0 (OU starts in its stationary law; band-limited is a random-phase sum) |
| Fractional (Kasdin) | **Q_d·H·Hᵀ**, with H[n,j] = h_{n−j} | the filter starts from rest, so not Toeplitz, even for α < 1 |

For the band-limited model, c is the covariance of *what the generator actually
produces*, Σ_k w_k cos(2πf_kτ), not the Ci formula. The expectation is then exact,
not merely a very good approximation.

## 3. The three expectations

**Autocovariance.** With y = x − m, where the mean m is removed with weight c
(0 for `:known`, 1 for `:trajectory`, 1/M for `:ensemble`):

    E[yₙ y_{n+k}] = Σ[n,n+k] − c·(Rₙ + R_{n+k})/N + c·G/N²,     Rₙ = Σ_m Σ[n,m],  G = Σ_n Rₙ

Average this over the N−k pairs. The middle term is why subtracting each record's
own mean biases the ACF downward: Rₙ/N is the covariance of xₙ with the record mean,
which is large when correlations are long (OU with T ~ τ_c, and above all 1/f noise).

**Periodogram.** E|Xₖ|² = Σ_{n,m} Σ[n,m] e^{−iωₖ(n−m)}.

- Toeplitz: this collapses to Σ_{|l|<N} (N−|l|) c(l) e^{−iωₖl}, the true spectrum
  smoothed by the Fejér kernel. That smoothing is **leakage**, and it is computed with one FFT.
- Fractional: using x = Hw gives E|Xₖ|² = Q_d Σ_{L=1}^{N} |Σ_{m<L} hₘ e^{−iωₖm}|². The
  truncated sums carry both the leakage and the start-from-rest transient.

**Variance profile.** Var[xₙ] = Σ[n,n]: constant for stationary models, and
Q_d Σ_{k≤n} h_k² for Kasdin, which grows ∝ ln n at α = 1.

The code never builds Σ. It uses row sums, diagonal averages and FFTs, costing
O(N·L) or O(N log N). **The tests do build Σ explicitly** for N = 13 and 14 and
check every fast formula against brute-force matrix algebra, for all models, all
three `demean` modes, and both record-length parities.

## 4. `validate_noise_model(model, ens)`

It runs four checks, each a z-score (measurement − exact expectation)/SE, where
every SE comes from the spread *across independent trajectories*:

| Check | Statistic | Why it is valid |
|-------|-----------|-----------------|
| mean at t = T | (x̄ − μ)/√(Var/M) | the M end values are i.i.d. |
| variance at t = T | (s² − Var)/(Var·√(2/(M−1))) | exact for Gaussian data |
| autocovariance | max over lags of \|z\| | per-lag SE from trajectory spread |
| periodogram | max over frequencies of \|z\| | per-bin SE from trajectory spread |

**Why "max |z|" needs a larger threshold.** Taking the worst of 300 lags, some
value will exceed 3 by chance. The threshold is Bonferroni-corrected,
z* = √2·erfc⁻¹(α/n), so the chance of *any* false alarm in a check is ≤ α
(default 10⁻³): 3.29 for one statistic, 4.65 for 300 lags.

**It has power.** The tests show that correct generators pass (OU, white, pink,
α = 0.5, band-limited), and that a wrong τ_c (0.7 instead of 0.5) or a 30% wrong
amplitude fails.

**The spectral slope** is reported but is not pass/fail: β̂ from the measured
periodogram next to the *same fit applied to the expected periodogram*. That
second number, not α, is the correct target for β̂.

## 5. How to use it

```julia
m   = FractionalNoiseModel(1.0, 1e-3)
ens = generate_ensemble(m, (0.0, 4.0), 1e-3, 1000; rng = Xoshiro(1))
r   = validate_noise_model(m, ens)       # prints a PASS/FAIL table
passed(r)                                # true/false for scripts and tests
r.acf.C, r.acf.expected                  # plot measured vs expected ACF
r.psd.S, r.psd.expected                  # plot measured vs expected spectrum
```

**Limitations.** The z-scores use normal approximations, so use M ≳ 100 trajectories.
Expectations exist for the four models in this package. A new model needs
`_stationary_acov` (if it is stationary) or its own `_acov_parts` and
`expected_periodogram` methods.