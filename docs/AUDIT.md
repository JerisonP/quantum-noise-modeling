# Audit of the previous code and notebook

This file records what was wrong in the pre-rewrite library and in
`Noise_Validation_revised_1.ipynb`, how each problem was confirmed, and where
the rewrite fixes it. Where a finding says "verified", the numbers come from
an independent NumPy re-implementation of the same generator and estimator
(`scripts/audit/` has the scripts).

## The one idea behind most of the bugs

> **Compare an estimator with the exact expectation of *that estimator* on
> *that generator*, not with an idealised quantity it only approximates.**

The old notebook compared finite-record, mean-subtracted, time-averaged
estimates with infinite-record stationary theory. Wherever that gap is
larger than the Monte Carlo error, a correct generator looks "wrong", or the
text quietly claims agreement that the numbers don't show. The rewrite computes
the exact expectation of each estimator (Step 3), so every comparison is
like-for-like.

## A. Results the notebook reports that contradict its own text

| # | Where | What the output shows | What the text claims | Cause | Verified |
|---|-------|-----------------------|----------------------|-------|----------|
| A1 | §4.3 `pink_acf`, §8 report | normalised-ACF RMSE **0.2365** (report: 0.227) | "typically below 0.05" | Library "theory" `Q_d Σ h_k h_{k+m}` is **not** the expectation of the estimator. The estimator averages over time a process that starts from rest, and it subtracts each trajectory's own mean. | RMSE vs exact expectation of the estimator: **0.0026** |
| A2 | §5.6 `ou_saveat_invariance` | C(τ_c) = 0.27–0.30 at every dt | "agree with the exact values within MC scatter" (theory 0.368) | Per-trajectory mean subtraction biases C(τ) down by ≈ 2τ_cσ²/T = **0.10** at T = 10 | demeaned 0.266, known-mean 0.361, theory 0.368 |
| A3 | §4.4 `pink_alpha_sweep` | 95% CI for α=0.8: [0.805, 0.833]; for α=1.5: [1.527, 1.558] | "every β̂ recovers α within its 95% CI" | Finite record started from rest, plus spectral leakage, steepens the expected periodogram | exact E[periodogram] slope: 0.8025 and 1.533, both inside the CIs |
| A4 | §7 `conv_pink` panel (c) | dashed target line drawn at α = 1 | table uses the binned-theory target 0.999 | plot and table use different targets | — |
| A5 | §4.2 `pink_running_stats` | 18× `ERROR: syntax error` | — | the GR backend cannot parse `\le` in `L"…"` labels | — |

## B. Numerical bugs in the library

| # | Function | Bug | Consequence | Fixed in |
|---|----------|-----|-------------|----------|
| B1 | `theoretical_autocorrelation(::FractionalNoiseModel)` for α<1 | Kasdin Eq. 110 evaluated with raw Γ functions | **NaN for every lag ≥ 172**. `validate_1f_noise.jl` uses 800 lags at α = 0.5 and 0.95, so those checks can't pass | Step 1 (stable recursion) |
| B2 | `noise_power_spectrum` (Wiener–Khinchin) | rectangle rule counts C(0) twice | constant offset 2Δτ·C(0); even with the *exact* OU ACF as input, output / theory = 3.5× at 5 Hz and 41× at 20 Hz | Step 2 |
| B3 | `noise_autocorrelation` | the default mean-subtraction mode is biased (A1, A2), and the method flipped between versions (the checkpoint used the grand mean, the current file per-trajectory means) | wrong ACF at long lags | Step 2 (explicit `mean=` option) + Step 3 |
| B4 | `load_propagator` | reads `[param, time, Λ…]`, but `save_propagator` writes `[time, Λ…]` | save → load round trip is misaligned or throws a BoundsError | Step 5 |
| B5 | `noise_autocorrelation` | triple loop O(M·N·L), about 5×10⁹ operations for the reference ensemble | slow | Step 2 (FFT, O(M N log N)) |

## C. Design problems (correct results, fragile code)

- **Three different PSD amplitude conventions.** The struct comment says `Q_psd/f^α`,
  `validate_1f_noise.jl` plots `Q_psd/(2πf)^α` (2× too low), and the code uses
  `2Q_psd/(2πf)^α`. The rewrite documents one: one-sided,
  S(f) ≈ 2Q_psd/(2πf)^α.
- **Irreproducible OU/white noise.** These ran through `EnsembleThreads` SDE
  solves seeded from the global RNG, so results depend on the thread count. White
  noise was faked with a stiff OU solve (τ_c = dt/10). The rewrite uses exact
  sampling with an explicit `rng` argument.
- **`dt` stored in `FractionalNoiseModel` as `Union{Float64,Nothing}`**, with
  warnings when it disagreed with the dt passed at run time. The rewrite passes dt
  only where it is used.
- Missing theory returned `NaN` (`TelegraphNoiseModel`, `PowerLawNoiseModel` in
  places), so a comparison against it silently turns into "N/A". The rewrite throws informative errors.
- `BandLimitedOneOverFNoiseModel` returns a `Vector` instead of the ensemble type.
- `Plots` is a hard dependency of the numerical core, and `_fit_psd_slope` is
  exported despite its underscore.
- The notebook's markdown numbers are stale: it says M = 400 (the code uses 1000), M_pool = 2000
  (6000), an R cap of 50 (60), R = 5 at M = 400 (25), and M = 100 in §4.5 (200), and it refers to "runtimes in
  Section 8" (they are in Section 9). `figures/` also contains files from earlier runs that no cell
  generates (`cmp_convergence.*`, `cmp_solution_convergence.*`).

## D. Checked and correct (kept)

- Kasdin pulse response h_k, AR coefficients a_k, driving variance
  Q_d = Q_psd·dt^(α−1), the discrete PSD 2Q_dΔt/|2sin(πfΔt)|^α with the one-sided factor
  of 2, and the finite-N variance Q_dΣh_k². All are re-derived and tested in Step 1.
- The OU SDE drift and diffusion, and the Lorentzian PSD.
- Physics: the gate phase φ(t) = ∫Ω/2, U₀ = exp(−iφσx), the column-major
  Liouvillian −i(I⊗H − Hᵀ⊗I), and average gate fidelity from the 6 Pauli eigenstates
  (they form a 2-design, so the 6-point average is exact).
