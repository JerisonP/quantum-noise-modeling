@testset "autocovariance: FFT result equals the defining sum (all demean modes)" begin
    rng = StableRNG(7)
    N, M, L = 37, 4, 10
    X = 0.3 .+ randn(rng, N, M)
    e = NoiseEnsemble((0:N-1) .* 0.1, X)
    for (mode, center) in ((:known, j -> 0.3), (:ensemble, j -> mean(X)), (:trajectory, j -> mean(X[:, j])))
        a = autocovariance(e; max_lag=L, demean=mode, μ=0.3)
        Cj = [sum((X[n, j] - center(j)) * (X[n+k, j] - center(j)) for n in 1:N-k) / (N - k)
              for k in 0:L, j in 1:M]
        @test a.C ≈ vec(mean(Cj; dims=2)) rtol = 1e-10
        @test a.C_err ≈ vec(std(Cj; dims=2)) ./ sqrt(M) rtol = 1e-8
        @test a.lags ≈ (0:L) .* 0.1
    end
    @test_throws ArgumentError autocovariance(e; max_lag=N)
    @test_throws ArgumentError autocovariance(e; demean=:median)
end

@testset "autocovariance: known/ensemble mean unbiased, per-trajectory mean biased" begin
    # OU with T = 10 = 20 τ_c. Per-trajectory demeaning subtracts ≈ Var(x̄_T) ≈ 2τ_cσ²/T = 0.1.
    m, dt = OUNoiseModel(1.0, 0.5), 0.01
    e = generate_ensemble(m, (0.0, 10.0), dt, 2000; rng=StableRNG(12))
    lagidx = (0, 25, 50, 100)
    for mode in (:known, :ensemble)
        a = autocovariance(e; max_lag=100, demean=mode)
        for k in lagidx
            @test within_se(a.C[k+1], theoretical_autocovariance(m, k * dt), a.C_err[k+1])
        end
    end
    a = autocovariance(e; max_lag=100, demean=:trajectory)
    k = 50                                                   # lag τ_c
    @test a.C[k+1] < theoretical_autocovariance(m, k * dt) - 5 * a.C_err[k+1]
    @test a.C[k+1] - theoretical_autocovariance(m, k * dt) ≈ -0.1 atol = 0.03
end

@testset "periodogram: definition, Parseval, white-noise level" begin
    rng = StableRNG(8)
    for N in (9, 10)                                         # odd and even (Nyquist bin present)
        X = randn(rng, N, 3)
        dt = 0.2
        p = periodogram(NoiseEnsemble((0:N-1) .* dt, X))
        K = N ÷ 2
        @test p.freqs ≈ (1:K) ./ (N * dt)
        brute = [2dt / N * abs2(sum(X[n, j] * cis(-2π * k * (n - 1) / N) for n in 1:N)) for k in 1:K, j in 1:3]
        @test p.S ≈ vec(mean(brute; dims=2)) rtol = 1e-10
        # Parseval: Σ w_k S_k df = mean over trajectories of the biased sample variance
        w = ones(K)
        iseven(N) && (w[K] = 0.5)
        @test sum(w .* p.S) / (N * dt) ≈ mean(var(X[:, j]; corrected=false) for j in 1:3) rtol = 1e-12
    end
    σ, dt = 1.5, 0.1
    p = periodogram(generate_ensemble(WhiteNoiseModel(σ), (0.0, 6.3), dt, 500; rng=StableRNG(9)))
    @test length(p.S) == 32
    # White-noise DFT bins are independent, each with sd ≈ its mean (√2× at Nyquist):
    # SE of the grand mean over 500 × 32 values ≤ √2 · 2σ²dt / √(500·32)
    @test within_se(mean(p.S), 2σ^2 * dt, sqrt(2) * 2σ^2 * dt / sqrt(500 * 32))
end

@testset "Wiener–Khinchin: exact for exact input; old formula off by 2Δτ·C(0)" begin
    m, dt = OUNoiseModel(1.3, 0.5), 0.01
    K = 5000                                                 # ϕ^K = e^{-100}: truncation negligible
    lags = (0:K) .* dt
    C = [theoretical_autocovariance(m, τ) for τ in lags]
    f = [0.0, 1.0, 7.3, 49.9]
    S = wiener_khinchin_psd(lags, C, f)
    @test S ≈ [theoretical_psd(m, x; dt) for x in f] rtol = 1e-9
    old = [4dt * sum(C .* cos.(2π .* x .* lags)) for x in f]  # previous implementation
    @test old .- S ≈ fill(2dt * C[1], length(f)) rtol = 1e-9
    r = wiener_khinchin_psd((lags=lags, C=C, C_err=zero(C)))
    @test length(r.freqs) == length(C) && r.freqs[end] ≈ 1 / (2dt)
    @test_throws ArgumentError wiener_khinchin_psd([0.1, 0.2], [1.0, 0.5], [1.0])   # lags must start at 0
end