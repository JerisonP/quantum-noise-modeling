# Where every formula comes from

Every formula in `src/` is listed here with its source. If a line of code
implements a formula, you should be able to find it in this table and then on
the cited page.

**Sources**
- **[K]** N. J. Kasdin, "Discrete simulation of colored noise and stochastic processes and 1/f^α power law noise generation", *Proc. IEEE* **83**, 802 (1995). §VI is on pp. 819–823.
- **[DH]** P. Dutta and P. M. Horn, "Low-frequency fluctuations in solids: 1/f noise", *Rev. Mod. Phys.* **53**, 497 (1981).
- **[RK]** J. Ruseckas and B. Kaulakys, "1/f noise from nonlinear stochastic differential equations", *Phys. Rev. E* **81**, 031105 (2010).
- **[W]** Wikipedia, "Ornstein–Uhlenbeck process" (sections *Definition*, *Formal solution*, *Numerical simulation*).

## Conventions

| Convention | Formula | Source |
|------------|---------|--------|
| One-sided PSD (package-wide) | S(f) = 4∫₀^∞ C(τ) cos(2πfτ) dτ, and C(0) = ∫₀^∞ S(f) df | [DH] p. 498 (unnumbered); [RK] Eq. 5 |
| Kasdin's PSD is two-sided | S_d(f) = Q_d Δt / (2 sin πfΔt)^α. White noise (α = 0) gives S = Q_dΔt over \|f\| ≤ 1/2Δt, which integrates to Q_d | [K] Eq. 98 |
| Hence | one-sided = 2 × Kasdin's S | — |

## OU noise (`src/noise/ou.jl`)

| Quantity | Formula | Source |
|----------|---------|--------|
| SDE | dζ = −(ζ−μ)/τ_c dt + σ√(2/τ_c) dW | [W] *Numerical simulation*; [W] *Definition* with θ = 1/τ_c, σ_W = σ√(2/τ_c) |
| Exact solution over one step | ζ(t+Δt) − μ = e^{−Δt/τ_c}(ζ(t) − μ) + (Gaussian, variance σ²(1−e^{−2Δt/τ_c})) | [W] *Formal solution* |
| Autocovariance | σ² e^{−\|τ\|/τ_c} | [W] stationary covariance |
| PSD | 4σ²τ_c / (1 + (2πfτ_c)²) | [DH] Eq. 7 (shape) + the one-sided convention above |

## 1/f^α noise (`src/noise/fractional.jl`)

| Quantity | Formula | Source |
|----------|---------|--------|
| Filter | H(z) = (1 − z⁻¹)^{−α/2} | [K] Eq. 97 (Hosking's fractional differencing) |
| α = 2 special case | H(z) = 1/(1 − z⁻¹), a random walk | [K] Eq. 96 |
| PSD (discrete, two-sided) | Q_dΔt / (2 sin πfΔt)^α | [K] Eq. 98 |
| PSD (low-f limit) | Q_d Δt^{1−α}/(2πf)^α | [K] Eq. 99 |
| Innovation variance | Q_d = Q / Δt^{1−α}, with `Q_psd` ≡ Q | [K] Eq. 100 |
| Pulse response (FIR) | h_0 = 1, h_k = (α/2 + k − 1) h_{k−1}/k | [K] Eqs. 102–104 |
| AR form | x_n = −a_1x_{n−1} − a_2x_{n−2} − … + w_n | [K] Eq. 115 |
| AR coefficients | a_0 = 1, a_k = (k − 1 − α/2) a_{k−1}/k | [K] Eq. 116 |
| Stationary variance (α < 1) | Q_d Γ(1−α)/Γ²(1−α/2) | [K] Eq. 111 |
| Stationary ACF (α < 1) | R_d(m) = Q_d Γ(α/2+m) / [2cos(πα/2) Γ(α) Γ(m+1−α/2)] = Q_d (−1)^m Γ(1−α)/[Γ(1+m−α/2) Γ(1−m−α/2)]. Code uses the ratio R_d(m)/R_d(m−1) = (m−1+α/2)/(m−α/2) | [K] Eqs. 109, 110 |
| Transient for α < 1 decays as t^{−α} | justifies "asymptotically stationary" | [K] footnote 6 |
| α = 1 nonstationary ACF | R(t,τ) ≅ (Q/2π)(ln 4t − ln\|τ\|) for t ≫ τ | [K] Eq. 93 |
| Finite-record variance | Var[x_n] = Q_d Σ_{k≤n} h_k² | follows from Eq. 104 with i.i.d. w |

## Band-limited 1/f (`src/noise/bandlimited.jl`)

| Quantity | Formula | Source |
|----------|---------|--------|
| PSD | Q/f on [f_l, f_h], 0 outside | your Mathematica model; the roll-offs are physically required ([DH] p. 499: ∫ f^{−α} df diverges at both ends) |
| ACF | Q[Ci(2πf_h\|τ\|) − Ci(2πf_l\|τ\|)] | ∫_{f_l}^{f_h} (Q/f) cos(2πfτ) df, the inverse of the one-sided convention above |
| Mathematica kernel | f_l = 1/2π, f_h = 100/2π gives Ci(100\|τ\|) − Ci(\|τ\|) | — |

Note: [DH] Eq. 8 with D(τ) ∝ 1/τ on [τ₁, τ₂] is the *soft* band-limited 1/f, a
sum of Lorentzians ([RK] Eq. 19 derives the same structure). Its ACF is
E₁(t/τ₂) − E₁(t/τ₁) rather than the Ci form. The brick-wall model here is the one
that matches your Mathematica kernel.

## Citations corrected from the old code

| Old citation | Actually |
|--------------|----------|
| Wiener–Khinchin "Dutta & Horn Eq. 5" | [DH] Eq. 5 is Voss's conditional mean. The WK formula is unnumbered on [DH] p. 498 and is [RK] Eq. 5 |
| "Dutta & Horn Eq. 93", R ∝ (Q/2π)(ln 4t − ln\|τ\|) | That is **[K] Eq. 93** |
| `PowerLawNoiseModel` weights "D(τ) ∝ τ^{α−1} from D&H Eq. 8" | [DH] give only D ∝ τ^{−1} (α = 1). The code's weights on a log-spaced τ grid amount to a variance density ∝ τ^{α−2}, which does give S ∝ f^{−α}, so the result was right but the citation was not |