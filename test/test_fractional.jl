const QNS = QuantumNoiseSimulator

@testset "constructor" begin
    m = FractionalNoiseModel(1, 1e-3)
    @test (m.α, m.Q_psd) === (1.0, 1e-3)
    @test_throws ArgumentError FractionalNoiseModel(-0.1, 1e-3)
    @test_throws ArgumentError FractionalNoiseModel(2.1, 1e-3)
    @test_throws ArgumentError FractionalNoiseModel(1.0, 0.0)
end

@testset "driving variance Q_d = Q_psd·dt^(α−1)" begin
    @test driving_variance(FractionalNoiseModel(1.0, 1e-3), 1e-2) ≈ 1e-3   # α = 1: dt-free
    @test driving_variance(FractionalNoiseModel(1.0, 1e-3), 1e-5) ≈ 1e-3
    @test driving_variance(FractionalNoiseModel(0.0, 1e-3), 1e-2) ≈ 1e-1   # white: Q/dt
    @test driving_variance(FractionalNoiseModel(2.0, 1e-3), 1e-2) ≈ 1e-5   # Brownian: Q·dt
    @test_throws ArgumentError driving_variance(FractionalNoiseModel(1.0, 1e-3), 0.0)
end

@testset "filter coefficients" begin
    @test pulse_response(1.0, 4) ≈ [1.0, 0.5, 0.375, 0.3125]      # h_k = h_{k−1}(α/2+k−1)/k
    @test pulse_response(0.0, 5) == [1.0, 0, 0, 0, 0]             # identity ⇒ white
    @test pulse_response(2.0, 5) == ones(5)                       # running sum ⇒ random walk
    @test ar_coefficients(2.0, 4) == [1.0, -1.0, 0.0, 0.0]        # (1 − z⁻¹)
    @test ar_coefficients(0.0, 4) == [1.0, 0.0, 0.0, 0.0]
    # h and a are inverse power series: (Σ h_k z⁻ᵏ)(Σ a_k z⁻ᵏ) = 1
    for α in (0.3, 1.0, 1.7)
        N = 40
        h, a = pulse_response(α, N), ar_coefficients(α, N)
        c = [sum(h[k+1] * a[n-k+1] for k in 0:n) for n in 0:N-1]
        @test c ≈ [1.0; zeros(N - 1)] atol = 1e-12
    end
    @test_throws ArgumentError pulse_response(1.0, 0)
end

@testset "fir_filter is linear causal convolution; FIR ≡ AR" begin
    rng = StableRNG(5)
    w = randn(rng, 50)
    h = pulse_response(1.3, 50)
    naive = [sum(h[k] * w[n-k+1] for k in 1:n) for n in 1:50]
    @test fir_filter(h, w) ≈ naive rtol = 1e-12
    @test fir_filter(pulse_response(2.0, 50), w) ≈ cumsum(w) rtol = 1e-12
    @test fir_filter(pulse_response(0.0, 50), w) ≈ w rtol = 1e-12
    W = randn(rng, 300, 3)
    for α in (0.5, 1.0, 1.5)
        @test fir_filter(pulse_response(α, 300), W) ≈ ar_filter(ar_coefficients(α, 300), W) rtol = 1e-9
    end
    @test_throws ArgumentError fir_filter(h[1:10], w)            # too few coefficients
end

@testset "generate_ensemble: shape, reproducibility, method equivalence" begin
    m = FractionalNoiseModel(1.0, 1e-3)
    e_fir = generate_ensemble(m, (0.0, 0.3), 1e-3, 3; rng=StableRNG(9))
    e_ar = generate_ensemble(m, (0.0, 0.3), 1e-3, 3; rng=StableRNG(9), method=:ar)
    @test size(samples(e_fir)) == (301, 3)
    @test samples(e_fir) ≈ samples(e_ar) rtol = 1e-9            # same innovations, two filters
    @test samples(generate_ensemble(m, (0.0, 0.3), 1e-3, 3; rng=StableRNG(9))) == samples(e_fir)
    @test_throws ArgumentError generate_ensemble(m, (0.0, 0.3), 1e-3, 3; method=:foo)
    # More trajectories than one internal batch (256) still works
    @test ntrajectories(generate_ensemble(m, (0.0, 0.01), 1e-3, 600; rng=StableRNG(1))) == 600
end

@testset "generated variance at t = T equals Q_d Σ h_k² (all α)" begin
    # The end-of-record cross-section is i.i.d. N(0, Q_d Σ_{k<N} h_k²), exactly.
    dt, N, M = 1e-3, 512, 4000
    for α in (0.0, 0.5, 1.0, 1.5, 2.0)
        m = FractionalNoiseModel(α, 1e-3)
        X = samples(generate_ensemble(m, (0.0, (N - 1) * dt), dt, M; rng=StableRNG(round(Int, 100α))))
        v = driving_variance(m, dt) * sum(abs2, pulse_response(α, N))
        @test within_se(mean(X[end, :]), 0.0, sqrt(v / M))
        @test within_se(var(X[end, :]), v, v * sqrt(2 / (M - 1)))
    end
end

@testset "stationary theory (α < 1)" begin
    m, dt = FractionalNoiseModel(0.5, 1e-3), 1e-3
    Qd = driving_variance(m, dt)
    # Recursion agrees with Kasdin Eq. 110 wherever Γ does not overflow
    for k in 0:50
        kasdin = Qd * (-1)^k * QNS.gamma(1 - m.α) /
                 (QNS.gamma(1 + k - m.α / 2) * QNS.gamma(1 - k - m.α / 2))
        @test theoretical_autocovariance(m, k * dt; dt) ≈ kasdin rtol = 1e-10
    end
    # ... and stays finite far beyond it (Eq. 110 gives NaN from k ≈ 172)
    @test isfinite(theoretical_autocovariance(m, 800dt; dt))
    @test theoretical_autocovariance(m, 0.0; dt) ≈ stationary_variance(m; dt)
    # Long-record limit of the FIR sum: Σ_k h_k h_{k+m} → R(m)/Q_d. The tail
    # h_k h_{k+m} ~ k^(α−2) makes the truncation error ~ N^(α−1) ≈ 1e-3 here,
    # hence the loose tolerance (this checks the formula, not the precision).
    h = pulse_response(m.α, 1_000_000)
    for lag in (0, 10)
        s = sum(h[k] * h[k+lag] for k in 1:length(h)-lag)
        @test Qd * s ≈ theoretical_autocovariance(m, lag * dt; dt) rtol = 5e-3
    end
    # ∫₀^{1/2dt} S(f) df = σ²  (substitute f = u² to remove the f^{-1/2} singularity)
    n = 200_000
    umax = sqrt(1 / (2dt))
    hu = umax / n
    us = ((1:n) .- 0.5) .* hu
    integral = sum(u -> 2u * theoretical_psd(m, u^2; dt), us) * hu
    @test integral ≈ stationary_variance(m; dt) rtol = 1e-5
    # White limit
    m0 = FractionalNoiseModel(0.0, 1e-3)
    @test stationary_variance(m0; dt) ≈ driving_variance(m0, dt)
    @test theoretical_autocovariance(m0, dt; dt) == 0.0
end

@testset "theory guards" begin
    pink = FractionalNoiseModel(1.0, 1e-3)
    @test_throws DomainError stationary_variance(pink; dt=1e-3)             # nonstationary
    @test_throws DomainError theoretical_autocovariance(pink, 0.0; dt=1e-3)
    @test_throws ArgumentError stationary_variance(FractionalNoiseModel(0.5, 1e-3))   # no dt
    @test_throws ArgumentError theoretical_autocovariance(FractionalNoiseModel(0.5, 1e-3), 1.5e-3; dt=1e-3)
end

@testset "PSD: discrete → continuous as f·dt → 0" begin
    for α in (0.5, 1.0, 1.5)
        m = FractionalNoiseModel(α, 1e-3)
        @test theoretical_psd(m, 1.0; dt=1e-6) ≈ theoretical_psd(m, 1.0) rtol = 1e-6
        @test theoretical_psd(m, 1.0) ≈ 2e-3 / (2π)^α
    end
    @test isinf(theoretical_psd(FractionalNoiseModel(1.0, 1e-3), 0.0))
end
