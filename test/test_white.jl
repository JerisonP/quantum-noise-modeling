@testset "constructor and theory" begin
    m = WhiteNoiseModel(0.5)
    @test (m.μ, m.σ) == (0.0, 0.5)
    @test_throws ArgumentError WhiteNoiseModel(-0.1)
    @test stationary_variance(m) == 0.25
    @test theoretical_autocovariance(m, 0.0) == 0.25
    @test theoretical_autocovariance(m, 0.01) == 0.0
    @test theoretical_autocovariance(m, 0.004; dt=0.01) == 0.25   # |τ| < dt/2 is lag 0
    @test theoretical_psd(m, 3.0; dt=0.01) ≈ 2 * 0.25 * 0.01
    @test_throws ArgumentError theoretical_psd(m, 3.0)            # needs dt
end

@testset "generated samples" begin
    μ, σ, M = -0.2, 1.5, 20_000
    X = samples(generate_ensemble(WhiteNoiseModel(μ, σ), (0.0, 1.0), 0.1, M; rng=StableRNG(3)))
    @test size(X) == (11, M)
    for k in (1, 6, 11)
        @test within_se(mean(X[k, :]), μ, σ / sqrt(M))
        @test within_se(var(X[k, :]), σ^2, σ^2 * sqrt(2 / (M - 1)))
    end
    prods = (X[3, :] .- μ) .* (X[4, :] .- μ)                      # independent in time
    @test within_se(mean(prods), 0.0, std(prods) / sqrt(M))
end