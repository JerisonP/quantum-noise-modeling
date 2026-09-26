# ─────────────────────────────────────────────────────────────────────────────
# The noise-model interface
# ─────────────────────────────────────────────────────────────────────────────
#
# Every noise model is a small immutable struct holding only *physical*
# parameters. The sampling interval dt is NOT stored in the model: it belongs to
# the simulation, and is passed whenever it matters. Theory functions follow
# one convention:
#
#     dt === nothing  →  the continuous-time process ζ(t)
#     dt given        →  exactly the sampled sequence ζ_k = ζ(k·dt) that
#                        generate_ensemble produces (aliasing included)
#
# Terminology: C(τ) = E[(ζ(t) − μ)(ζ(t+τ) − μ)] is the AUTOCOVARIANCE (it has
# units of ζ²). The normalised ρ(τ) = C(τ)/C(0) is the autocorrelation.

"""
    AbstractNoiseModel

Supertype of all noise models. A concrete model must implement

- `generate_ensemble(model, tspan, dt, n_samples; rng)`, which returns a [`NoiseEnsemble`](@ref)

and, where they exist in closed form,

- `noise_mean(model)`, which defaults to `0.0`
- `stationary_variance(model; dt=nothing)`
- `theoretical_autocovariance(model, τ; dt=nothing)`
- `theoretical_psd(model, f; dt=nothing)`, the **one-sided** PSD, with ∫₀^∞ S(f) df = variance

A method that has no closed form throws an informative error rather than
returning `NaN`, so a missing theory can never silently pass a check.
"""
abstract type AbstractNoiseModel end

"""
    generate_ensemble(model, tspan, dt, n_samples; rng=Random.default_rng()) -> NoiseEnsemble

Draw `n_samples` independent trajectories of `model` on the uniform grid
`tspan[1] : dt : tspan[2]`. The grid must contain an integer number of steps
of size `dt`. Pass an explicit `rng` (e.g. `Xoshiro(1)`) for reproducibility.
"""
function generate_ensemble end

"Mean E[ζ(t)] of the process."
noise_mean(::AbstractNoiseModel) = 0.0

"""
    stationary_variance(model; dt=nothing)

Stationary variance σ² = C(0). Throws a `DomainError` for nonstationary
models, where the variance depends on the record length.
"""
function stationary_variance end

"""
    theoretical_autocovariance(model, τ; dt=nothing)

Stationary autocovariance C(τ) = E[(ζ(t)−μ)(ζ(t+τ)−μ)] at lag `τ`.
"""
function theoretical_autocovariance end

"""
    theoretical_psd(model, f; dt=nothing)

One-sided power spectral density S(f), with f in cycles per unit time.
When `dt` is given, this is the exact PSD of the sampled sequence on
0 ≤ f ≤ 1/(2dt) (and it integrates to the variance over that band).
"""
function theoretical_psd end

# ── Shared argument checking ─────────────────────────────────────────────────

"""
    _time_grid(tspan, dt) -> Vector{Float64}

`tspan[1], tspan[1]+dt, …, tspan[2]`. Errors if the span is not an integer
number of steps, because silently rounding would change the record length T.
"""
function _time_grid(tspan, dt::Real)
    t0, t1 = Float64(tspan[1]), Float64(tspan[2])
    dt > 0 || throw(ArgumentError("dt must be > 0, got $dt"))
    t1 > t0 || throw(ArgumentError("tspan must satisfy tspan[2] > tspan[1], got $tspan"))
    steps = (t1 - t0) / dt
    n = round(Int, steps)
    isapprox(steps, n; atol=1e-6) || throw(ArgumentError(
        "(tspan[2] - tspan[1]) / dt = $steps is not an integer; choose dt dividing the span"))
    return collect(t0 .+ (0:n) .* Float64(dt))
end

function _check_nsamples(n_samples::Integer)
    n_samples >= 1 || throw(ArgumentError("n_samples must be ≥ 1, got $n_samples"))
    return nothing
end
