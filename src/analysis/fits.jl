# ─────────────────────────────────────────────────────────────────────────────
# Regression utilities: turn a measured curve into a parameter ± error
# ─────────────────────────────────────────────────────────────────────────────

"""
    student_t_quantile(p, ν) -> t

Quantile of Student's t distribution with ν degrees of freedom, via the
inverse regularised incomplete beta function: for p > 1/2,
P(|T| > t) = I_{ν/(ν+t²)}(ν/2, 1/2) = 2(1−p).
"""
function student_t_quantile(p::Real, ν::Real)
    0 < p < 1 || throw(ArgumentError("p must lie in (0, 1), got $p"))
    ν > 0 || throw(ArgumentError("ν must be > 0, got $ν"))
    p == 0.5 && return 0.0
    p < 0.5 && return -student_t_quantile(1 - p, ν)
    x = first(beta_inc_inv(ν / 2, 0.5, 2 * (1 - p)))
    return sqrt(ν * (1 - x) / x)
end

"""
    ols(x, y; level=0.95) -> (intercept, slope, se_intercept, se_slope, slope_ci, r2, n)

Ordinary least squares y = a + b·x, with the standard errors from the residual
variance s² = RSS/(n−2) and a t-based `level` confidence interval for the slope.
The errors are valid when the residuals are independent with a common variance.
"""
function ols(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; level::Real=0.95)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("x and y differ in length"))
    n >= 3 || throw(ArgumentError("need ≥ 3 points for standard errors, got $n"))
    x̄, ȳ = mean(x), mean(y)
    Sxx = sum(abs2, x .- x̄)
    Sxx > 0 || throw(ArgumentError("x has no spread"))
    b = sum((x .- x̄) .* (y .- ȳ)) / Sxx
    a = ȳ - b * x̄
    rss = sum(abs2, y .- (a .+ b .* x))
    tss = sum(abs2, y .- ȳ)
    s2 = rss / (n - 2)
    se_b = sqrt(s2 / Sxx)
    se_a = sqrt(s2 * (1 / n + x̄^2 / Sxx))
    tq = student_t_quantile((1 + level) / 2, n - 2)
    return (intercept=a, slope=b, se_intercept=se_a, se_slope=se_b,
            slope_ci=(b - tq * se_b, b + tq * se_b),
            r2=tss > 0 ? 1 - rss / tss : 1.0, n=n)
end

"""
    powerlaw_fit(freqs, S; fmin=nothing, fmax=nothing, margin=0.3, nbins=0, level=0.95)
        -> (β, se, ci, log10_amplitude, r2, n, fmin, fmax)

Fit S(f) = A·f^(−β) by OLS of log₁₀S on log₁₀f.

**Fit band.** Give `fmin`/`fmax` explicitly, or else the central part of the
log-frequency range is used, dropping a fraction `margin` of the decades at
each end (the old `_fit_psd_slope` convention).

**`nbins > 0`: log-binned fit.** Average S inside `nbins` equal-width bins in
log f, then fit the bin means. On a linear frequency grid, most raw points sit
in the top decade, so an unbinned fit is dominated by it. Binning gives every
decade equal weight and reduces the scatter of β̂.

`se`/`ci` are regression errors, valid only if residuals are independent. For
a statement about the *estimator*, use replicate ensembles (Step 3/notebook
§7), and compare β̂ with the same fit applied to the exact expected spectrum,
not with α.
"""
function powerlaw_fit(freqs::AbstractVector{<:Real}, S::AbstractVector{<:Real};
                      fmin=nothing, fmax=nothing, margin::Real=0.3, nbins::Integer=0,
                      level::Real=0.95)
    length(freqs) == length(S) || throw(DimensionMismatch("freqs and S differ in length"))
    ok = (freqs .> 0) .& (S .> 0)
    lf, lS = log10.(freqs[ok]), log10.(S[ok])
    lo = fmin === nothing ? lf[1] + margin * (lf[end] - lf[1]) : log10(fmin)
    hi = fmax === nothing ? lf[end] - margin * (lf[end] - lf[1]) : log10(fmax)
    band = (lf .>= lo) .& (lf .<= hi)
    fband = freqs[ok][band]                  # the frequencies actually used (reported as-is)
    lf, lS = lf[band], lS[band]
    if nbins > 0
        edges = range(minimum(lf), maximum(lf); length=nbins + 1)
        xs, ys = Float64[], Float64[]
        for b in 1:nbins
            inb = (lf .>= edges[b]) .& (b == nbins ? (lf .<= edges[b+1]) : (lf .< edges[b+1]))
            if any(inb)
                push!(xs, mean(lf[inb]))
                push!(ys, log10(mean(10 .^ lS[inb])))   # average S in the bin, then log
            end
        end
        lf, lS = xs, ys
    end
    length(lf) >= 3 || throw(ArgumentError("fewer than 3 points in the fit band"))
    fit = ols(lf, lS; level)
    return (β=-fit.slope, se=fit.se_slope, ci=(-fit.slope_ci[2], -fit.slope_ci[1]),
            log10_amplitude=fit.intercept, r2=fit.r2, n=fit.n,
            fmin=minimum(fband), fmax=maximum(fband))
end

"""
    exponential_fit(lags, C; threshold=0.05, level=0.95) -> (τ_c, se, ci, C0, r2, n)

Fit C(τ) = C₀·e^{−τ/τ_c} by OLS of ln C on τ, using the leading run of lags
with C > `threshold`·C(0). The run stops at the first lag that fails, so noisy
tail values that happen to be positive are never included. The standard error
of τ_c = −1/slope comes from the delta method: se(τ_c) = se(slope)/slope².
"""
function exponential_fit(lags::AbstractVector{<:Real}, C::AbstractVector{<:Real};
                         threshold::Real=0.05, level::Real=0.95)
    length(lags) == length(C) || throw(DimensionMismatch("lags and C differ in length"))
    C[1] > 0 || throw(ArgumentError("C(0) must be positive"))
    n = something(findfirst(c -> !(c > threshold * C[1]), C), length(C) + 1) - 1
    n >= 3 || throw(ArgumentError("fewer than 3 lags above threshold"))
    fit = ols(lags[1:n], log.(C[1:n]); level)
    b = fit.slope
    b < 0 || throw(ArgumentError("autocovariance is not decaying (slope = $b)"))
    τ = -1 / b
    se = fit.se_slope / b^2
    tq = student_t_quantile((1 + level) / 2, n - 2)
    return (τ_c=τ, se=se, ci=(τ - tq * se, τ + tq * se), C0=exp(fit.intercept), r2=fit.r2, n=n)
end