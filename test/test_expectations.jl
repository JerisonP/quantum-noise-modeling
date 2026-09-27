const QNSe = QuantumNoiseSimulator

# ── Brute force: build Σ explicitly and apply the estimator definitions ─────
"Explicit covariance matrix Σ[n,m] = Cov(x_n, x_m) of the sampled generator output."
function brute_sigma(m, N, dt)
    if m isa FractionalNoiseModel
        h = pulse_response(m.α, N)
        H = zeros(N, N)
        for j in 1:N, n in j:N
            H[n, j] = h[n-j+1]
        end
        return driving_variance(m, dt) .* (H * H')
    end
    c = QNSe._stationary_acov(m, N, dt)
    return [c[abs(n - k)+1] for n in 1:N, k in 1:N]
end

"E[Ĉ(k)] from Σ: y = x − c·(sample mean), averaged over the N − k pairs."
function brute_acov(Σ, L, c)
    N = size(Σ, 1)
    R = vec(sum(Σ; dims=2))
    Σy = Σ .- c .* (R .+ R') ./ N .+ c * sum(R) / N^2
    return [mean(Σy[n, n+k] for n in 1:N-k) for k in 0:L]
end

"E[Ŝ(f_k)] = (2dt/N) Σ_{n,m} Σ[n,m] e^{−iω_k(n−m)}."
brute_pgram(Σ, dt) = (N = size(Σ, 1);
    [2dt / N * real(sum(cis(-2π * k * (n - m) / N) * Σ[n, m] for n in 1:N, m in 1:N)) for k in 1:N÷2])

@testset "fast formulas equal brute-force matrix algebra" begin
    dt = 0.02
    models = (OUNoiseModel(0.4, 1.7, 0.35), WhiteNoiseModel(0.8),
              BandLimitedOneOverFNoiseModel(1.0, 1 / (2π), 100 / (2π), 64),
              FractionalNoiseModel(0.5, 1e-3), FractionalNoiseModel(1.0, 1e-3),
              FractionalNoiseModel(1.5, 1e-3))
    for m in models, N in (13, 14)                     # odd and even record lengths
        Σ = brute_sigma(m, N, dt)
        scale = maximum(abs, Σ)
        for (mode, c) in ((:known, 0.0), (:trajectory, 1.0), (:ensemble, 1 / 7))
            E = expected_autocovariance(m, N, dt; max_lag=N - 1, demean=mode, M=7)
            @test E.C ≈ brute_acov(Σ, N - 1, c) atol = 1e-12 * scale
        end
        @test expected_periodogram(m, N, dt).S ≈ brute_pgram(Σ, dt) atol = 1e-12 * scale * N * dt
        @test expected_variance(m, N, dt) ≈ [Σ[n, n] for n in 1:N]
    end
end

@testset "consistency with Step 1 theory" begin
    # Stationary + known mean: the expectation is just the theoretical autocovariance
    m = OUNoiseModel(1.3, 0.5)
    E = expected_autocovariance(m, 1001, 0.01; max_lag=200)
    @test E.C ≈ [theoretical_autocovariance(m, τ) for τ in E.lags]
    @test expected_variance(m, 50, 0.01) == fill(1.3^2, 50)
    # Fractional: variance profile Q_d Σ h², and lag-0 known-mean = its time average
    f = FractionalNoiseModel(1.0, 1e-3)
    v = expected_variance(f, 512, 1e-3)
    @test v[end] ≈ driving_variance(f, 1e-3) * sum(abs2, pulse_response(1.0, 512))
    @test issorted(v)                                     # nonstationary: grows with time
    @test expected_autocovariance(f, 512, 1e-3; max_lag=0).C[1] ≈ mean(v)
    # Argument checking
    @test_throws ArgumentError expected_autocovariance(m, 100, 0.01; demean=:ensemble)  # needs M
    @test_throws ArgumentError expected_autocovariance(m, 100, 0.01; demean=:median)
    @test_throws ArgumentError expected_autocovariance(m, 100, 0.01; max_lag=100)
end

@testset "old notebook §4.3: pink-noise ACF agrees once compared with the right expectation" begin
    # α = 1, per-trajectory demeaning (the old estimator). Against the old "theory"
    # Q_d Σ h_k h_{k+m} the normalised ACF was off by RMSE ≈ 0.24; against the exact
    # expectation of the same estimator it agrees to Monte Carlo precision.
    m, dt, N, L = FractionalNoiseModel(1.0, 1e-3), 1e-3, 1001, 300
    e = generate_ensemble(m, (0.0, (N - 1) * dt), dt, 2000; rng=StableRNG(31))
    a = autocovariance(e; max_lag=L, demean=:trajectory)
    E = expected_autocovariance(m, e; max_lag=L, demean=:trajectory)
    h = pulse_response(1.0, N)
    old = [sum(h[i] * h[i+k] for i in 1:N-k) for k in 0:L]
    ρ = a.C ./ a.C[1]
    rmse(x) = sqrt(mean(abs2, ρ .- x ./ x[1]))
    @test rmse(E.C) < 0.02                                # ≈ 0.004
    @test rmse(old) > 0.1                                 # ≈ 0.27: the old failure, reproduced
end

@testset "old notebook §4.4: the right target for the fitted spectral slope" begin
    # The fit of the EXACT expected periodogram (N = 4001, dt = 1e-3, log-binned,
    # 24 bins, 30% margins) is the value β̂ should scatter around. It differs from α
    # because of finite-record leakage and the start-from-rest transient.
    for (α, target) in ((0.8, 0.802499), (1.5, 1.533201))
        Ep = expected_periodogram(FractionalNoiseModel(α, 1e-3), 4001, 1e-3)
        @test powerlaw_fit(Ep.freqs, Ep.S; nbins=24).β ≈ target atol = 2e-4
    end
end