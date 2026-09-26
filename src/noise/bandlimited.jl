# ─────────────────────────────────────────────────────────────────────────────
# Band-limited stationary 1/f noise (matches the Mathematica kernel)
# ─────────────────────────────────────────────────────────────────────────────

"""
    BandLimitedOneOverFNoiseModel(Q_psd, fl, fh, extension_factor=8)

Stationary Gaussian noise with a one-sided PSD that is exactly 1/f inside a band
and zero outside it:

    S(f) = Q_psd / f   for fl ≤ f ≤ fh,     0 otherwise      (f in cycles / unit time)

# Exact statistics
- variance        C(0) = Q_psd · ln(fh/fl)
- autocovariance  C(τ) = Q_psd · [Ci(2π fh |τ|) − Ci(2π fl |τ|)],  where Ci is the cosine integral

**Matching the Mathematica kernel.** Choosing `fl = 1/(2π)` and `fh = 100/(2π)`
gives C(τ) = Q_psd·[Ci(100|τ|) − Ci(|τ|)]:

```julia
m = BandLimitedOneOverFNoiseModel(Q, 1/(2π), 100/(2π))
```

Unlike [`FractionalNoiseModel`](@ref) with α = 1, this process is **stationary**,
because the band cutoff at `fl` removes the divergent low-frequency power.

# How it is generated
Random-phase Fourier synthesis on a periodic trace `L ≈ extension_factor × n_t`
samples long, done with one inverse FFT per trajectory. The first `n_t` samples are kept:

    x(t) = Σ_k √w_k · (a_k cos 2πf_k t + b_k sin 2πf_k t),   a_k, b_k ~ N(0, 1)

- f_k = k/(L·dt) is the FFT frequency grid, and w_k = ∫ S(f) df over bin k (the exact
  bin integral, not S(f_k)·df). This makes the variance exactly Q_psd·ln(fh/fl).
- The synthesised process is exactly Gaussian and stationary, with
  C_gen(τ) = Σ_k w_k cos(2πf_k τ). That is a quadrature of the Ci formula, accurate
  to ≲ 10⁻³·C(0) for the default `extension_factor = 8`.
- **Why the extension:** the synthesis is periodic with period L·dt, so C_gen(τ)
  is only right for τ ≪ L·dt. With `extension_factor = 1` the ACF would be wrong
  by O(1) at lags comparable to the record.

Errors are thrown if `fh` is at or above Nyquist (the band cannot be represented)
or if the frequency resolution 1/(L·dt) exceeds `fl` (the bottom of the band cannot be resolved).
"""
struct BandLimitedOneOverFNoiseModel <: AbstractNoiseModel
    Q_psd::Float64
    fl::Float64
    fh::Float64
    extension_factor::Int

    function BandLimitedOneOverFNoiseModel(Q_psd::Real, fl::Real, fh::Real,
                                           extension_factor::Integer=8)
        Q_psd > 0 || throw(ArgumentError("Q_psd must be > 0, got $Q_psd"))
        fl > 0 || throw(ArgumentError("fl must be > 0, got $fl"))
        fh > fl || throw(ArgumentError("fh must be > fl, got fl = $fl, fh = $fh"))
        extension_factor >= 1 || throw(ArgumentError("extension_factor must be ≥ 1"))
        return new(Float64(Q_psd), Float64(fl), Float64(fh), Int(extension_factor))
    end
end

stationary_variance(m::BandLimitedOneOverFNoiseModel; dt=nothing) = m.Q_psd * log(m.fh / m.fl)

function theoretical_autocovariance(m::BandLimitedOneOverFNoiseModel, τ::Real; dt=nothing)
    x = abs(τ)
    iszero(x) && return stationary_variance(m)
    return m.Q_psd * (cosint(2π * m.fh * x) - cosint(2π * m.fl * x))
end

"""
    theoretical_psd(m::BandLimitedOneOverFNoiseModel, f; dt=nothing)

Q_psd/f inside [fl, fh] and 0 outside. The band lies below Nyquist (this is enforced
at generation), so the sampled sequence has the same PSD and `dt` is ignored.
"""
function theoretical_psd(m::BandLimitedOneOverFNoiseModel, f::Real; dt=nothing)
    return m.fl <= f <= m.fh ? m.Q_psd / f : 0.0
end

"""
    _synthesis_grid(m, n_t, dt) -> (f, w, k, L)

The frequency bins used by `generate_ensemble`: FFT indices `k`, frequencies
`f = k/(L·dt)`, and weights `w` = exact ∫ S(f) df over each bin ∩ [fl, fh].
`sum(w) == Q_psd·ln(fh/fl)`. The generated process has autocovariance
`Σ w .* cos.(2π f τ)` exactly; the validation layer uses that.
"""
function _synthesis_grid(m::BandLimitedOneOverFNoiseModel, n_t::Integer, dt::Real)
    L = nextprod([2, 3, 5], m.extension_factor * n_t)
    df = 1 / (L * dt)
    kcap = cld(L, 2) - 1                       # highest bin that is not DC or Nyquist
    m.fh <= (kcap + 0.5) * df || throw(ArgumentError(
        "fh = $(m.fh) is at or above the Nyquist frequency 1/(2dt) = $(1 / (2dt)); reduce dt"))
    df <= m.fl || throw(ArgumentError(
        "frequency resolution 1/(L·dt) = $df exceeds fl = $(m.fl); increase " *
        "extension_factor (now $(m.extension_factor)) or the record length"))
    ks = Int[]
    ws = Float64[]
    for k in max(1, floor(Int, m.fl / df - 0.5)):min(kcap, ceil(Int, m.fh / df + 0.5))
        lo = max((k - 0.5) * df, m.fl)
        hi = min((k + 0.5) * df, m.fh)
        if hi > lo
            push!(ks, k)
            push!(ws, m.Q_psd * log(hi / lo))
        end
    end
    return (f=ks .* df, w=ws, k=ks, L=L)
end

function generate_ensemble(m::BandLimitedOneOverFNoiseModel, tspan, dt::Real, n_samples::Integer;
                           rng::AbstractRNG=default_rng())
    t = _time_grid(tspan, dt)
    _check_nsamples(n_samples)
    n_t = length(t)
    grid = _synthesis_grid(m, n_t, dt)
    L = grid.L
    # irfft(Y, L)[n] = (1/L)[Y₀ + 2 Re Σ_k Y_k e^{2πikn/L}], so Y_k = (L/2)·√w_k·(a_k − i b_k)
    # gives x_n = Σ_k √w_k (a_k cos + b_k sin), as documented.
    amp = sqrt.(grid.w) .* (L / 2)
    Y = zeros(ComplexF64, L ÷ 2 + 1)
    P = plan_irfft(Y, L)
    X = Matrix{Float64}(undef, n_t, n_samples)
    for j in 1:n_samples
        fill!(Y, 0)
        for (i, k) in enumerate(grid.k)
            a = randn(rng)
            b = randn(rng)
            Y[k+1] = amp[i] * complex(a, -b)
        end
        x = P * Y
        X[:, j] .= view(x, 1:n_t)
    end
    return NoiseEnsemble(t, X)
end