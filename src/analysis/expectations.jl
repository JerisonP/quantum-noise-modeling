# ─────────────────────────────────────────────────────────────────────────────
# Exact finite-record expectations of the estimators
# ─────────────────────────────────────────────────────────────────────────────
#
# This file is the "answer key" for Step 2's estimators. For a record of n_t
# samples at spacing dt, each function returns exactly what the estimator
# would give on average over infinitely many trajectories:
#
#   E[autocovariance(ens; demean)]   → expected_autocovariance
#   E[periodogram(ens)]              → expected_periodogram
#   Var[x_n] at every time point     → expected_variance
#
# Everything follows from the covariance matrix Σ[n,m] = Cov(x_n, x_m) of the
# sampled generator output:
#
#   • OU, white, band-limited: stationary, so Σ is Toeplitz, Σ[n,m] = c(|n−m|),
#     where c is the exact covariance of what the generator produces.
#   • FractionalNoiseModel: the filter starts from rest, so x = H·w with H lower
#     triangular (H[n,j] = h_{n−j}) and Σ = Q_d·H·Hᵀ. This is NOT Toeplitz, even
#     for α < 1, and for α ≥ 1 it is strongly nonstationary.
#
# The old notebook compared estimates with infinite-record stationary formulas;
# comparing with these exact expectations is what fixes its three contradicted results.
#
# Derivation of the demeaned autocovariance. With y = x − m and
#   m = μ (known), the trajectory's own mean (c = 1), or the grand mean over M trajectories (c = 1/M):
#     E[y_n y_{n+k}] = Σ[n,n+k] − c·(R_n + R_{n+k})/N + c·G/N²,
#   where R_n = Σ_m Σ[n,m] (row sums) and G = Σ_n R_n. Averaging over the N−k
#   pairs gives E[Ĉ(k)]. The test suite checks every formula here against brute-force
#   matrix algebra.

# ── Covariance structure per model ───────────────────────────────────────────

"""
    _stationary_acov(model, n_t, dt) -> c

`c[l+1]` = exact Cov(x_n, x_{n+l}) of the sampled generator output, for
l = 0 … n_t−1. It is defined for stationary models only.
"""
_stationary_acov(m::AbstractNoiseModel, n_t, dt) = throw(ArgumentError(
    "no finite-record expectation is implemented for $(typeof(m))"))

_stationary_acov(m::OUNoiseModel, n_t, dt) =
    [theoretical_autocovariance(m, l * dt) for l in 0:n_t-1]

_stationary_acov(m::WhiteNoiseModel, n_t, dt) = [l == 0 ? m.σ^2 : 0.0 for l in 0:n_t-1]

function _stationary_acov(m::BandLimitedOneOverFNoiseModel, n_t, dt)
    g = _synthesis_grid(m, n_t, dt)          # exactly the bins the generator uses
    return [sum(g.w .* cos.(2π .* g.f .* (l * dt))) for l in 0:n_t-1]
end

"""
    _acov_parts(model, N, dt, L) -> (diag, R)

`diag[k+1]` = (1/(N−k)) Σ_n Σ[n,n+k] for k = 0…L (the known-mean expectation), and
`R[n+1]` = Σ_m Σ[n,m] (row sums, needed for mean subtraction).
"""
function _acov_parts(m::AbstractNoiseModel, N, dt, L)
    c = _stationary_acov(m, N, dt)
    cs = cumsum(c)
    R = [cs[n+1] + cs[N-n] - c[1] for n in 0:N-1]   # Σ_{l=0}^{n} c_l + Σ_{l=1}^{N−1−n} c_l
    return c[1:L+1], R
end

function _acov_parts(m::FractionalNoiseModel, N, dt, L)
    h = pulse_response(m.α, N)
    Qd = driving_variance(m, dt)
    # Σ[n,n+k] = Q_d Σ_{i=0}^{n} h_i h_{i+k}; averaging over n = 0…N−1−k gives
    # (Q_d/(N−k)) Σ_i (N−k−i) h_i h_{i+k}
    diag = zeros(L + 1)
    for k in 0:L
        A = 0.0
        B = 0.0
        @inbounds for i in 0:N-1-k
            p = h[i+1] * h[i+k+1]
            A += p
            B += i * p
        end
        diag[k+1] = Qd * ((N - k) * A - B) / (N - k)
    end
    # Row sums R = Q_d·H·(Hᵀ·1); (Hᵀ·1)_j = Σ_{m=0}^{N−1−j} h_m is a reversed cumulative sum
    R = Qd .* fir_filter(h, reverse(cumsum(h)))
    return diag, R
end

# ── Public API ───────────────────────────────────────────────────────────────

"""
    expected_autocovariance(model, n_t, dt; max_lag, demean=:known, M=nothing) -> (lags, C)
    expected_autocovariance(model, ens; max_lag, demean=:known)

Exact expectation of [`autocovariance`](@ref) for records of `n_t` samples at
spacing `dt`, with the same `demean` convention:

- `:known`: E = mean of Σ's k-th diagonal. For a stationary model this is simply C(k·dt).
- `:trajectory`: includes the downward bias from subtracting each record's own mean.
- `:ensemble`: includes the (small) bias of the grand mean, which needs `M`.

The ensemble form reads `n_t`, `dt` and `M` from `ens`.
"""
function expected_autocovariance(m::AbstractNoiseModel, n_t::Integer, dt::Real;
                                 max_lag::Integer=floor(Int, 0.3 * (n_t - 1)),
                                 demean::Symbol=:known, M=nothing)
    N = Int(n_t)
    0 <= max_lag <= N - 1 || throw(ArgumentError("max_lag must lie in 0:$(N-1), got $max_lag"))
    c = if demean === :known
        0.0
    elseif demean === :trajectory
        1.0
    elseif demean === :ensemble
        M === nothing && throw(ArgumentError("demean = :ensemble needs the number of trajectories M"))
        1.0 / M
    else
        throw(ArgumentError("demean must be :known, :ensemble or :trajectory, got :$demean"))
    end
    diag, R = _acov_parts(m, N, Float64(dt), max_lag)
    G = sum(R)
    cR = [0.0; cumsum(R)]                            # cR[i+1] = Σ_{j<i} R_j
    C = similar(diag)
    for k in 0:max_lag
        n = N - k
        s = (cR[n+1] - cR[1]) + (cR[N+1] - cR[k+1])  # Σ_{i<n} R_i + Σ_{i≥k} R_i
        C[k+1] = diag[k+1] - c * s / (N * n) + c * G / N^2
    end
    return (lags=collect(0:max_lag) .* dt, C=C)
end

expected_autocovariance(m::AbstractNoiseModel, e::NoiseEnsemble; kwargs...) =
    expected_autocovariance(m, ntimes(e), timestep(e); M=ntrajectories(e), kwargs...)

"""
    expected_periodogram(model, n_t, dt) -> (freqs, S)
    expected_periodogram(model, ens)

Exact expectation of [`periodogram`](@ref), including finite-record leakage and,
for `FractionalNoiseModel`, the start-from-rest transient. Comparing a measured
spectrum (or a fitted slope) with this, rather than with the ideal S(f), is
what makes the α-sweep consistent.

- Stationary models: E|X_k|² = Σ_{|l|<N} (N−|l|) c(l) e^{−iω_k l} (the Fejér-kernel
  smoothing of the true spectrum), evaluated with one FFT.
- `FractionalNoiseModel`: E|X_k|² = Q_d Σ_{L=1}^{N} |Σ_{m<L} h_m e^{−iω_k m}|².
"""
function expected_periodogram(m::AbstractNoiseModel, n_t::Integer, dt::Real)
    N = Int(n_t)
    c = _stationary_acov(m, N, dt)
    # Fold negative lags onto the length-N DFT: lag l−N ≡ l (mod N)
    a = [l == 0 ? N * c[1] : (N - l) * c[l+1] + l * c[N-l+1] for l in 0:N-1]
    F = rfft(a)
    K = N ÷ 2
    return (freqs=collect(1:K) ./ (N * dt), S=(2dt / N) .* real.(F[2:K+1]))
end

function expected_periodogram(m::FractionalNoiseModel, n_t::Integer, dt::Real)
    N = Int(n_t)
    h = pulse_response(m.α, N)
    Qd = driving_variance(m, dt)
    K = N ÷ 2
    S = zeros(K)
    for k in 1:K
        ω = 2π * k / N
        G = zero(ComplexF64)
        acc = 0.0
        @inbounds for n in 0:N-1
            G += h[n+1] * cis(-ω * n)                # partial DFT sum G_{n+1}(ω)
            acc += abs2(G)
        end
        S[k] = (2dt / N) * Qd * acc
    end
    return (freqs=collect(1:K) ./ (N * dt), S=S)
end

expected_periodogram(m::AbstractNoiseModel, e::NoiseEnsemble) =
    expected_periodogram(m, ntimes(e), timestep(e))

"""
    expected_variance(model, n_t, dt) -> Vector

Var[x_n] at each of the `n_t` time points: constant for stationary models, and
Q_d Σ_{k≤n} h_k² (growing) for `FractionalNoiseModel`.
"""
expected_variance(m::AbstractNoiseModel, n_t::Integer, dt::Real) =
    fill(first(_stationary_acov(m, n_t, dt)), n_t)

expected_variance(m::FractionalNoiseModel, n_t::Integer, dt::Real) =
    driving_variance(m, dt) .* cumsum(abs2.(pulse_response(m.α, n_t)))