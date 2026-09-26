# ─────────────────────────────────────────────────────────────────────────────
# Ornstein–Uhlenbeck noise
# ─────────────────────────────────────────────────────────────────────────────

"""
    OUNoiseModel(μ, σ, τ_c)
    OUNoiseModel(σ, τ_c)        # zero-mean shorthand

Ornstein–Uhlenbeck process, the standard model of noise with a single
correlation time τ_c:

    dζ = −(ζ − μ)/τ_c dt + σ √(2/τ_c) dW

This is the parameterisation in Wikipedia's *Ornstein–Uhlenbeck process →
Numerical simulation* section (stationary std σ, correlation time τ_c). The
*Definition* section's form dx = θ(μ − x)dt + σ_W dW maps to it through

    θ = 1/τ_c,     σ_W = σ √(2/τ_c)     (stationary variance σ_W²/(2θ) = σ²).

# Exact statistics (stationary)
- mean μ, variance σ²
- autocovariance  C(τ) = σ² exp(−|τ|/τ_c)
- one-sided PSD   S(f) = 4σ²τ_c / (1 + (2πfτ_c)²), a Lorentzian with corner f_c = 1/(2πτ_c).
  This is S(f) = 4∫₀^∞ C(τ) cos(2πfτ) dτ (Dutta & Horn p. 498, Ruseckas & Kaulakys
  Eq. 5), and the shape is Dutta & Horn Eq. 7.

# How it is generated (exactly, with no integrator error)
Wikipedia's *Formal solution* section solves the SDE exactly:
x_t = μ + (x_s − μ)e^{−θ(t−s)} + σ_W ∫_s^t e^{−θ(t−u)} dW_u. The stochastic
integral is Gaussian with mean 0 and variance σ²(1 − e^{−2(t−s)/τ_c}), and it is
independent of x_s. Applying that solution over each grid step gives the AR(1)
recursion

    ζ_{k+1} − μ = ϕ (ζ_k − μ) + σ √(1 − ϕ²) ξ_k,     ϕ = exp(−Δt/τ_c),  ξ_k ~ N(0, 1)

with ζ_0 ~ N(μ, σ²) drawn from the stationary law. This is not an
approximation: the sampled sequence has *exactly* the OU finite-dimensional
distributions for any Δt. There is therefore no time-step error to study,
which is why the old "convergence in dt" experiment is replaced by an exactness test.
"""
struct OUNoiseModel <: AbstractNoiseModel
    μ::Float64
    σ::Float64
    τ_c::Float64

    function OUNoiseModel(μ::Real, σ::Real, τ_c::Real)
        σ >= 0 || throw(ArgumentError("σ must be ≥ 0, got $σ"))
        τ_c > 0 || throw(ArgumentError("τ_c must be > 0, got $τ_c"))
        return new(Float64(μ), Float64(σ), Float64(τ_c))
    end
end

OUNoiseModel(σ::Real, τ_c::Real) = OUNoiseModel(0.0, σ, τ_c)

function generate_ensemble(m::OUNoiseModel, tspan, dt::Real, n_samples::Integer;
                           rng::AbstractRNG=default_rng())
    t = _time_grid(tspan, dt)
    _check_nsamples(n_samples)
    n_t = length(t)
    ϕ = exp(-dt / m.τ_c)
    s = m.σ * sqrt(-expm1(-2dt / m.τ_c))    # σ√(1−ϕ²); expm1 keeps precision when dt ≪ τ_c
    X = Matrix{Float64}(undef, n_t, n_samples)
    @inbounds for j in 1:n_samples
        x = m.σ * randn(rng)                 # stationary initial condition
        X[1, j] = m.μ + x
        for k in 2:n_t
            x = ϕ * x + s * randn(rng)
            X[k, j] = m.μ + x
        end
    end
    return NoiseEnsemble(t, X)
end

noise_mean(m::OUNoiseModel) = m.μ

stationary_variance(m::OUNoiseModel; dt=nothing) = m.σ^2

theoretical_autocovariance(m::OUNoiseModel, τ::Real; dt=nothing) =
    m.σ^2 * exp(-abs(τ) / m.τ_c)

"""
    theoretical_psd(m::OUNoiseModel, f; dt=nothing)

- `dt === nothing`: continuous-time Lorentzian 4σ²τ_c / (1 + (2πfτ_c)²).
- `dt` given: exact one-sided PSD of the sampled AR(1) sequence,
  2Δt σ²(1−ϕ²) / (1 − 2ϕ cos(2πfΔt) + ϕ²). This is the Lorentzian with
  aliasing included; the two agree when f ≪ 1/(2Δt).
"""
function theoretical_psd(m::OUNoiseModel, f::Real; dt=nothing)
    if dt === nothing
        return 4 * m.σ^2 * m.τ_c / (1 + (2π * f * m.τ_c)^2)
    end
    ϕ = exp(-dt / m.τ_c)
    return 2 * dt * m.σ^2 * (1 - ϕ^2) / (1 - 2ϕ * cos(2π * f * dt) + ϕ^2)
end