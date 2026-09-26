# ─────────────────────────────────────────────────────────────────────────────
# White (i.i.d. Gaussian) noise
# ─────────────────────────────────────────────────────────────────────────────

"""
    WhiteNoiseModel(μ, σ)
    WhiteNoiseModel(σ)          # zero-mean shorthand

Independent Gaussian samples ζ_k ~ N(μ, σ²) on the simulation grid.

**Note on units.** `σ` is the per-sample standard deviation, so the physical
noise strength depends on Δt: the one-sided PSD of the sequence is
S = 2σ²Δt (flat up to Nyquist). To hold the *PSD* fixed while changing Δt,
scale σ ∝ 1/√Δt. Continuous white noise has no finite variance, so
`theoretical_psd` requires `dt`.

(The previous version faked white noise by running an OU SDE solver with
τ_c = Δt/10. That is slow and stiff, and it is only approximately white; direct
sampling is exact.)
"""
struct WhiteNoiseModel <: AbstractNoiseModel
    μ::Float64
    σ::Float64

    function WhiteNoiseModel(μ::Real, σ::Real)
        σ >= 0 || throw(ArgumentError("σ must be ≥ 0, got $σ"))
        return new(Float64(μ), Float64(σ))
    end
end

WhiteNoiseModel(σ::Real) = WhiteNoiseModel(0.0, σ)

function generate_ensemble(m::WhiteNoiseModel, tspan, dt::Real, n_samples::Integer;
                           rng::AbstractRNG=default_rng())
    t = _time_grid(tspan, dt)
    _check_nsamples(n_samples)
    X = randn(rng, length(t), n_samples)
    X .= m.μ .+ m.σ .* X
    return NoiseEnsemble(t, X)
end

noise_mean(m::WhiteNoiseModel) = m.μ

stationary_variance(m::WhiteNoiseModel; dt=nothing) = m.σ^2

"""
    theoretical_autocovariance(m::WhiteNoiseModel, τ; dt=nothing)

σ² at zero lag, 0 otherwise. With `dt` given, any |τ| < dt/2 counts as lag 0.
"""
function theoretical_autocovariance(m::WhiteNoiseModel, τ::Real; dt=nothing)
    is_zero_lag = dt === nothing ? iszero(τ) : abs(τ) < dt / 2
    return is_zero_lag ? m.σ^2 : 0.0
end

function theoretical_psd(m::WhiteNoiseModel, f::Real; dt=nothing)
    dt === nothing && throw(ArgumentError(
        "sampled white noise has no continuous-time PSD; pass dt (S = 2σ²dt)"))
    return 2 * m.σ^2 * dt
end
