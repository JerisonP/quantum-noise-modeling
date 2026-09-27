# ─────────────────────────────────────────────────────────────────────────────
# validate_noise_model: measurement vs exact expectation, with honest error bars
# ─────────────────────────────────────────────────────────────────────────────

"""
    ValidationCheck(name, statistic, threshold, passed, note)

One pass/fail check. `statistic` is a |z|-score, (estimate − expectation)/SE,
or its maximum over many lags or frequencies. It passes when it is at most `threshold`.
"""
struct ValidationCheck
    name::String
    statistic::Float64
    threshold::Float64
    passed::Bool
    note::String
end

"""
    ValidationReport

Returned by [`validate_noise_model`](@ref). Fields:

- `model`, `n_t`, `M`, `dt`: what was validated
- `checks::Vector{ValidationCheck}`: the pass/fail checks
- `acf = (lags, C, C_err, expected)`: measured autocovariance vs its exact expectation
- `psd = (freqs, S, S_err, expected)`: measured periodogram vs its exact expectation
- `slope = (β, β_expected, se)`: power-law fits to both spectra (informational)

`passed(report)` is `true` when every check passes. The data fields are there so
the notebook can plot exactly what was tested.
"""
struct ValidationReport{Mo<:AbstractNoiseModel}
    model::Mo
    n_t::Int
    M::Int
    dt::Float64
    checks::Vector{ValidationCheck}
    acf::NamedTuple
    psd::NamedTuple
    slope::NamedTuple
end

"`true` if every check in the report passed."
passed(r::ValidationReport) = all(c -> c.passed, r.checks)

"""
    _bonferroni_z(n, α)

Threshold z* with P(max of n independent |N(0,1)| > z*) ≈ α, i.e. a family-wise
false-alarm rate α across n simultaneous checks. Examples: n = 1, α = 1e-3 gives 3.29;
n = 300 gives 4.65.
"""
_bonferroni_z(n, α) = sqrt(2) * erfcinv(α / n)

"max |(a − b)/se| over entries, treating 0/0 as agreement."
function _max_abs_z(a, b, se)
    z = 0.0
    for i in eachindex(a, b, se)
        d = abs(a[i] - b[i])
        d == 0 && continue
        z = max(z, se[i] > 0 ? d / se[i] : Inf)
    end
    return z
end

"""
    validate_noise_model(model, ens; max_lag, demean=:known, α=1e-3, nbins=24, margin=0.3)
        -> ValidationReport

Test whether the ensemble `ens` has the statistics of `model`. Every
comparison is against the **exact finite-record expectation of the same
estimator** (see `expected_*`), so a correct generator passes at any record length or α,
and the checks do not depend on idealised stationary formulas.

Checks, each with family-wise false-alarm rate `α`:

1. **mean at t = T**: the M end-of-record values are i.i.d., so z = (x̄ − μ)/√(Var/M).
2. **variance at t = T**: z = (s² − Var)/(Var·√(2/(M−1))). This is exact for Gaussian noise.
3. **autocovariance**: the maximum over lags of |Ĉ(k) − E Ĉ(k)|/SE(k), where SE comes from the
   spread across trajectories.
4. **periodogram**: the maximum over frequencies of |Ŝ(f) − E Ŝ(f)|/SE(f).

Not a pass/fail check: `slope` compares the power-law fit of the measured
periodogram with the same fit applied to the expected periodogram (the right
target for β̂, not α itself).

The z-scores use normal approximations, so use M ≳ 100 trajectories.
"""
function validate_noise_model(model::AbstractNoiseModel, ens::NoiseEnsemble;
                              max_lag::Integer=floor(Int, 0.3 * (ntimes(ens) - 1)),
                              demean::Symbol=:known, α::Real=1e-3,
                              nbins::Integer=24, margin::Real=0.3)
    X = samples(ens)
    N, M = size(X)
    M >= 2 || throw(ArgumentError("need at least 2 trajectories"))
    dt = timestep(ens)
    μ = noise_mean(model)
    checks = ValidationCheck[]

    # 1–2: end-of-record cross-section (i.i.d. across trajectories)
    v = expected_variance(model, N, dt)[end]
    xT = X[end, :]
    z1 = _bonferroni_z(1, α)
    zm = abs(mean(xT) - μ) / sqrt(v / M)
    push!(checks, ValidationCheck("mean at t = T", zm, z1, zm <= z1,
                                  "x̄ = $(round(mean(xT); sigdigits=4)), expected $(round(μ; sigdigits=4))"))
    s2 = sum(abs2, xT .- mean(xT)) / (M - 1)
    zv = abs(s2 - v) / (v * sqrt(2 / (M - 1)))
    push!(checks, ValidationCheck("variance at t = T", zv, z1, zv <= z1,
                                  "s² = $(round(s2; sigdigits=4)), expected $(round(v; sigdigits=4))"))

    # 3: autocovariance vs its exact expectation
    a = autocovariance(ens; max_lag, demean, μ)
    Ea = expected_autocovariance(model, N, dt; max_lag, demean, M)
    za = _max_abs_z(a.C, Ea.C, a.C_err)
    ta = _bonferroni_z(max_lag + 1, α)
    push!(checks, ValidationCheck("autocovariance ($(max_lag + 1) lags, demean = :$demean)", za, ta,
                                  za <= ta, "max |z| over lags"))

    # 4: periodogram vs its exact expectation
    p = periodogram(ens)
    Ep = expected_periodogram(model, N, dt)
    zp = _max_abs_z(p.S, Ep.S, p.S_err)
    tp = _bonferroni_z(length(p.S), α)
    push!(checks, ValidationCheck("periodogram ($(length(p.S)) frequencies)", zp, tp, zp <= tp,
                                  "max |z| over frequencies"))

    # Informational: spectral slope vs the same fit on the expected spectrum
    slope = try
        fm = powerlaw_fit(p.freqs, p.S; nbins, margin)
        fe = powerlaw_fit(Ep.freqs, Ep.S; nbins, margin)
        (β=fm.β, β_expected=fe.β, se=fm.se)
    catch
        (β=NaN, β_expected=NaN, se=NaN)
    end

    return ValidationReport(model, N, M, dt, checks,
                            (lags=a.lags, C=a.C, C_err=a.C_err, expected=Ea.C),
                            (freqs=p.freqs, S=p.S, S_err=p.S_err, expected=Ep.S),
                            slope)
end

function Base.show(io::IO, ::MIME"text/plain", r::ValidationReport)
    println(io, "Noise validation: ", r.model)
    println(io, "  ", r.M, " trajectories × ", r.n_t, " points, dt = ", r.dt)
    for c in r.checks
        mark = c.passed ? "PASS" : "FAIL"
        println(io, "  [", mark, "] ", rpad(c.name, 44), " |z| = ",
                rpad(round(c.statistic; digits=2), 6), " (limit ", round(c.threshold; digits=2), ")  ",
                c.note)
    end
    s = r.slope
    println(io, "  spectral slope: β̂ = ", round(s.β; digits=4), " ± ", round(s.se; digits=4),
            ", expected ", round(s.β_expected; digits=4))
    print(io, "  overall: ", passed(r) ? "PASS" : "FAIL")
end

Base.show(io::IO, r::ValidationReport) =
    print(io, "ValidationReport(", typeof(r.model).name.name, ", ", passed(r) ? "PASS" : "FAIL", ")")