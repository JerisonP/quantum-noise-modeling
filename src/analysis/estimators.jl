# ─────────────────────────────────────────────────────────────────────────────
# Estimators: what we MEASURE from an ensemble
# ─────────────────────────────────────────────────────────────────────────────
#
# Every estimator returns its per-trajectory spread as a standard error. The
# trajectories are independent, so SE = std(per-trajectory estimates)/√M is
# an honest error bar, even though samples WITHIN a trajectory are correlated.
#
# Each estimator is defined exactly, so that Step 3 can compute its exact
# expectation for every noise model. The notebook then compares measurements
# with that expectation, never with an idealised quantity the estimator only
# approximates.

"""
    autocovariance(ens; max_lag, demean=:known, μ=0.0) -> (lags, C, C_err)

Time-averaged autocovariance, averaged over trajectories:

    Ĉ(k) = (1/M) Σ_j  (1/(N−k)) Σ_{n=1}^{N−k} y_{n,j} y_{n+k,j},     y = x − m

for lags k = 0 … `max_lag` (default 30% of the record). This is the "unbiased"
normalisation 1/(N−k). The choice of `m` matters, and it is explicit:

| `demean`      | m                          | when to use |
|---------------|----------------------------|-------------|
| `:known`      | the true mean `μ`          | the mean is known (all our models) → **no bias from mean estimation** |
| `:ensemble`   | grand mean over all samples | unknown mean; its bias ~ 1/(M·T) is negligible for large M |
| `:trajectory` | each trajectory's own mean  | the old behaviour. It is biased low by ≈ Var(x̄_T), which is ≈ 2τ_cσ²/T for OU and much worse for 1/f noise |

`C_err` is the standard error across trajectories. Computed by zero-padded
FFT in O(M·N log N).
"""
function autocovariance(e::NoiseEnsemble; max_lag::Integer=floor(Int, 0.3 * (ntimes(e) - 1)),
                        demean::Symbol=:known, μ::Real=0.0)
    X = samples(e)
    N, M = size(X)
    0 <= max_lag <= N - 1 || throw(ArgumentError("max_lag must lie in 0:$(N-1), got $max_lag"))
    grand = demean === :known ? Float64(μ) :
            demean === :ensemble ? mean(X) :
            demean === :trajectory ? NaN :
            throw(ArgumentError("demean must be :known, :ensemble or :trajectory, got :$demean"))
    nfft = nextprod([2, 3, 5], 2N - 1)             # zero padding ⇒ linear (not circular) correlation
    buf = zeros(nfft)
    Pf = plan_rfft(buf)
    Pi = plan_irfft(Pf * buf, nfft)
    Cj = Matrix{Float64}(undef, max_lag + 1, M)    # per-trajectory estimates
    inv_pairs = [1.0 / (N - k) for k in 0:max_lag]   # 1/(N−k) normalisation
    for j in 1:M
        x = view(X, :, j)
        m = demean === :trajectory ? mean(x) : grand
        fill!(buf, 0.0)
        buf[1:N] .= x .- m
        F = Pf * buf
        r = Pi * complex.(abs2.(F))                # r[k+1] = Σ_n y_n y_{n+k}
        Cj[:, j] .= view(r, 1:max_lag+1) .* inv_pairs
    end
    C = vec(mean(Cj; dims=2))
    C_err = M > 1 ? vec(std(Cj; dims=2)) ./ sqrt(M) : fill(NaN, max_lag + 1)
    return (lags=collect(0:max_lag) .* timestep(e), C=C, C_err=C_err)
end

"""
    periodogram(ens) -> (freqs, S, S_err)

Ensemble-averaged **one-sided** periodogram at the Fourier frequencies
f_k = k/(N·dt), for k = 1 … ⌊N/2⌋ (DC is excluded):

    Ŝ(f_k) = (2 dt / N) · |Σ_n x_n e^{−2πikn/N}|²

Every bin, the Nyquist bin included, is a PSD *density*: for white noise
E[Ŝ] = 2σ²dt in every bin. Normalisation check (Parseval, tested): with bin
widths df = 1/(N dt) and half weight on the Nyquist bin when N is even,
Σ_k w_k Ŝ(f_k) df equals each trajectory's (biased) sample variance, averaged
over trajectories. Subtracting a constant from x changes only the excluded DC
bin, so no demeaning is needed. No window is applied: the leakage this causes
is part of the estimator, and Step 3 accounts for it exactly.
"""
function periodogram(e::NoiseEnsemble)
    X = samples(e)
    N, M = size(X)
    dt = timestep(e)
    K = N ÷ 2
    K >= 1 || throw(ArgumentError("need at least 2 time points"))
    c = 2dt / N
    P = plan_rfft(zeros(N))
    Sj = Matrix{Float64}(undef, K, M)
    for j in 1:M
        F = P * X[:, j]
        @inbounds for k in 1:K
            Sj[k, j] = c * abs2(F[k+1])
        end
    end
    S = vec(mean(Sj; dims=2))
    S_err = M > 1 ? vec(std(Sj; dims=2)) ./ sqrt(M) : fill(NaN, K)
    return (freqs=collect(1:K) ./ (N * dt), S=S, S_err=S_err)
end

"""
    wiener_khinchin_psd(lags, C, freqs) -> S
    wiener_khinchin_psd(acf; freqs=range(0, 1/(2Δτ); length=length(acf.C)))

One-sided PSD from an autocovariance, S(f) = 4∫₀^∞ C(τ) cos(2πfτ) dτ
(Dutta & Horn p. 498; Ruseckas & Kaulakys Eq. 5), discretised by the
**trapezoid rule** on the lag grid τ_k = kΔτ:

    S(f) = 2Δτ [ C_0 + 2 Σ_{k≥1} C_k cos(2πfτ_k) ]

For the exact autocovariance of a sampled sequence (and enough lags) this
equals that sequence's exact PSD. The old implementation used 4Δτ Σ_{k≥0},
which counts C_0 twice and adds a constant 2Δτ·C_0 at every frequency (AUDIT B2).
Truncating at a finite maximum lag causes ripple, and a truncated
*estimated* ACF is biased at f ≲ 1/τ_max.
"""
function wiener_khinchin_psd(lags::AbstractVector{<:Real}, C::AbstractVector{<:Real},
                             freqs::AbstractVector{<:Real})
    length(lags) == length(C) || throw(DimensionMismatch("lags and C differ in length"))
    length(lags) >= 2 || throw(ArgumentError("need at least 2 lags"))
    iszero(lags[1]) || throw(ArgumentError("lags must start at 0"))
    Δτ = lags[2] - lags[1]
    all(k -> isapprox(lags[k+1] - lags[k], Δτ; rtol=1e-6), 1:length(lags)-1) ||
        throw(ArgumentError("lags must be uniformly spaced"))
    return [2Δτ * (C[1] + 2 * sum(C[k] * cos(2π * f * lags[k]) for k in 2:length(C))) for f in freqs]
end

function wiener_khinchin_psd(acf::NamedTuple;
                             freqs=range(0, 1 / (2 * (acf.lags[2] - acf.lags[1])); length=length(acf.C)))
    return (freqs=collect(freqs), S=wiener_khinchin_psd(acf.lags, acf.C, freqs))
end