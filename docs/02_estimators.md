# Step 2: Estimators and fits

*How we measure statistics from an ensemble, and why each choice matters.*

## 1. The one rule

An estimator is a **recipe applied to finite data**, and its average over many
ensembles (its *expectation*) is generally **not** the idealised quantity it is
named after. Three effects cause the difference, and each broke a result in the old notebook:

| Effect | What it does | Old symptom |
|--------|--------------|-------------|
| Mean subtraction | Subtracting the sample mean removes Var(x̄), which is real variance | OU C(τ_c) read 0.27 instead of 0.37 (old §5.6) |
| Finite record, start from rest | 1/f^α noise has not "forgotten" t = 0 | 1/f ACF RMSE 0.24 (old §4.3) |
| Leakage (no window) | Power from low frequencies spills into high ones | α-sweep CIs missed α (old §4.4) |

So every estimator here is **defined exactly**, and Step 3 computes its exact
expectation. This step builds the measuring instruments, and Step 3 builds the
answer key for each one.

## 2. Autocovariance: `autocovariance(ens; max_lag, demean, μ)`

    Ĉ(k) = (1/M) Σ_j  (1/(N−k)) Σ_n (x_{n,j} − m)(x_{n+k,j} − m)

* **1/(N−k)**: there are N−k pairs at lag k, so this is a true average over pairs.
* **`demean` is explicit**, because it is the single biggest source of bias:

  - `:known` (default) subtracts the true μ, so there is no bias. All our models have a known mean.
  - `:ensemble` subtracts the grand mean over all M·N samples. The bias ≈ Var(grand mean) ~ 1/(M·T), negligible for large M.
  - `:trajectory` subtracts each path's own mean (the old behaviour). The bias ≈ −Var(x̄_T):
    ≈ −2τ_cσ²/T for OU, and much larger for 1/f noise, where the variance of the time average barely shrinks with T.

  The test suite *demonstrates* this. On OU with T = 20τ_c, `:known` and `:ensemble`
  sit within 1–2 SE of σ²e^{−τ/τ_c}, while `:trajectory` is low by 0.100 ≈ 20 SE.
  That is exactly the old §5.6 table.

* **Error bars.** Samples within one trajectory are correlated, so a naive
  SE from all N·M products would be far too small. Trajectories *are*
  independent, so we compute Ĉ separately for each one and take
  SE = std/√M. The same approach is used for every estimator in the package.

* **Speed.** The inner sum for all lags at once is an FFT correlation
  (zero-padded to avoid wrap-around), which is O(N log N) instead of O(N·L).

## 3. Periodogram: `periodogram(ens)`

    Ŝ(f_k) = (2Δt/N) |DFT(x)_k|²,     f_k = k/(NΔt),  k = 1…⌊N/2⌋

* **Factor 2Δt/N.** Δt/N turns |DFT|² into a density per unit frequency, and
  the 2 folds negative frequencies onto positive ones (the one-sided convention).
* **Parseval (tested exactly).** Σ Ŝ_k·df equals the sample variance, with half weight
  on the Nyquist bin, which sits at the edge of the band.
* **No demeaning needed.** A constant only affects the DC bin, which is dropped.
* **No window, on purpose.** For steep spectra (α ≥ 1) leakage lifts the high
  frequencies. Rather than hide this with a taper, Step 3 computes E[Ŝ]
  *including* leakage, which is why the α = 1.5 fit target becomes 1.533.

## 4. Wiener–Khinchin: `wiener_khinchin_psd`

The papers define the one-sided PSD as S(f) = 4∫₀^∞ C(τ) cos(2πfτ) dτ
([DH] p. 498, [RK] Eq. 5). The **trapezoid rule** on lags kΔτ gives

    S(f) = 2Δτ [C₀ + 2 Σ_{k≥1} C_k cos(2πfτ_k)]

where the endpoint C₀ gets half weight. The old code gave C₀ full weight (4Δτ Σ_{k≥0}),
which adds 2Δτ·C₀ at *every* frequency. That is harmless where S is large and
dominant where S is small: 41× too big at 20 Hz for OU. With the exact OU
covariance as input, the corrected formula reproduces the exact sampled-OU
spectrum to 10⁻⁹ (tested).

## 5. Fits

| Function | Model | Notes |
|----------|-------|-------|
| `ols(x, y)` | y = a + bx | SEs from s² = RSS/(n−2) and a t-based CI; checked against scipy's `linregress` |
| `student_t_quantile(p, ν)` | — | via the inverse incomplete beta function, so no Distributions.jl dependency; checked against scipy |
| `powerlaw_fit(f, S; nbins)` | S = A f^−β | explicit band or margin; `nbins > 0` gives equal weight per decade |
| `exponential_fit(lags, C)` | C₀e^{−τ/τ_c} | uses only the *leading run* above threshold, so noisy tail points never enter; delta-method SE |

**A caution on regression error bars.** The SE from a fit assumes independent
residuals. Periodogram ordinates are roughly independent, but ACF values at
neighbouring lags are strongly correlated, so an `exponential_fit` SE is too
optimistic. The honest error bar on τ̂_c comes from replicate ensembles (old
notebook §7, rebuilt in Step 6). The fit SE is reported but not relied on.