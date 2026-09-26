@testset "constructor" begin
    m = OUNoiseModel(0.3, 2.0)
    @test (m.μ, m.σ, m.τ_c) == (0.0, 0.3, 2.0)
    @test OUNoiseModel(1, 2, 3).μ === 1.0                # Int arguments are promoted
    @test_throws ArgumentError OUNoiseModel(-1.0, 1.0)
    @test_throws ArgumentError OUNoiseModel(1.0, 0.0)
end

@testset "closed-form theory" begin
    m = OUNoiseModel(1.5, 2.0, 0.5)
    @test noise_mean(m) == 1.5
    @test stationary_variance(m) == 4.0
    @test theoretical_autocovariance(m, 0.0) == 4.0
    @test theoretical_autocovariance(m, 0.5) ≈ 4.0 / ℯ
    @test theoretical_autocovariance(m, -0.3) == theoretical_autocovariance(m, 0.3)
    @test theoretical_psd(m, 0.0) ≈ 4 * 4.0 * 0.5                       # white plateau 4σ²τ_c
    @test theoretical_psd(m, 1 / (2π * 0.5)) ≈ theoretical_psd(m, 0.0) / 2   # corner frequency
end

@testset "conventions: Wikipedia parameters and one-sided Wiener–Khinchin" begin
    m = OUNoiseModel(1.3, 0.5)
    θ, σ_W = 1 / m.τ_c, m.σ * sqrt(2 / m.τ_c)          # Wikipedia *Definition* parameters
    @test σ_W^2 / (2θ) ≈ stationary_variance(m)         # Wikipedia stationary variance
    # S(f) = 4∫₀^∞ C(τ)cos(2πfτ)dτ  (Dutta & Horn p. 498; Ruseckas & Kaulakys Eq. 5)
    τ = range(0, 40 * m.τ_c; length=400_001)
    for f in (0.0, 0.3, 2.0)
        S = 4 * trapz(τ, [theoretical_autocovariance(m, x) * cos(2π * f * x) for x in τ])
        @test S ≈ theoretical_psd(m, f) rtol = 1e-6
    end
end

@testset "sampled PSD: integrates to σ² and matches the Lorentzian at low f" begin
    m, dt = OUNoiseModel(1.3, 0.5), 0.01
    f = range(0, 1 / (2dt); length=200_001)
    S = [theoretical_psd(m, x; dt) for x in f]
    @test trapz(f, S) ≈ 1.3^2 rtol = 1e-6                 # Parseval: exact variance
    @test theoretical_psd(m, 0.1; dt=1e-5) ≈ theoretical_psd(m, 0.1) rtol = 1e-6
end

@testset "generated samples: shape, grid, reproducibility" begin
    m = OUNoiseModel(1.0, 0.5)
    e1 = generate_ensemble(m, (0.0, 2.0), 0.01, 7; rng=StableRNG(1))
    e2 = generate_ensemble(m, (0.0, 2.0), 0.01, 7; rng=StableRNG(1))
    e3 = generate_ensemble(m, (0.0, 2.0), 0.01, 7; rng=StableRNG(2))
    @test size(samples(e1)) == (201, 7)
    @test times(e1) ≈ 0.0:0.01:2.0
    @test samples(e1) == samples(e2)
    @test samples(e1) != samples(e3)
    @test_throws ArgumentError generate_ensemble(m, (0.0, 1.0), 0.01, 0)
end

@testset "generated samples: exact stationary statistics" begin
    # Cross-sections X[k, :] are i.i.d. over trajectories, so their standard
    # errors are exact: SE(mean) = σ/√M, SE(var) ≈ σ²√(2/(M−1)).
    μ, σ, τ_c, dt, M = 0.7, 2.0, 0.5, 0.01, 20_000
    m = OUNoiseModel(μ, σ, τ_c)
    X = samples(generate_ensemble(m, (0.0, 3.0), dt, M; rng=StableRNG(11)))
    for k in (1, 151, 301)                               # t = 0, 1.5, 3: stationary from t = 0
        @test within_se(mean(X[k, :]), μ, σ / sqrt(M))
        @test within_se(var(X[k, :]), σ^2, σ^2 * sqrt(2 / (M - 1)))
    end
    # Two-time covariance at lag τ_c: E[(x_k−μ)(x_{k+m}−μ)] = σ²/e
    lag = round(Int, τ_c / dt)
    prods = (X[51, :] .- μ) .* (X[51+lag, :] .- μ)
    @test within_se(mean(prods), σ^2 / ℯ, std(prods) / sqrt(M))
    # Exactness of the AR(1) step: lag-1 regression coefficient equals ϕ
    y0, y1 = X[100, :] .- μ, X[101, :] .- μ
    ϕ̂ = sum(y0 .* y1) / sum(y0 .^ 2)
    @test within_se(ϕ̂, exp(-dt / τ_c), sqrt((1 - exp(-2dt / τ_c)) / M))
end