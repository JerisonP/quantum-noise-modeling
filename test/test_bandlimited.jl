const QNSb = QuantumNoiseSimulator
# Mathematica kernel: C(τ) = Q[Ci(100|τ|) − Ci(|τ|)]
mathematica_model(Q=1.0; ext=8) = BandLimitedOneOverFNoiseModel(Q, 1 / (2π), 100 / (2π), ext)

"Exact autocovariance of what generate_ensemble produces: Σ_k w_k cos(2π f_k τ)."
generated_acov(g, τ) = sum(g.w .* cos.(2π .* g.f .* τ))

@testset "constructor" begin
    m = mathematica_model()
    @test m.extension_factor == 8
    @test_throws ArgumentError BandLimitedOneOverFNoiseModel(0.0, 0.1, 1.0)
    @test_throws ArgumentError BandLimitedOneOverFNoiseModel(1.0, 0.0, 1.0)
    @test_throws ArgumentError BandLimitedOneOverFNoiseModel(1.0, 1.0, 0.5)
    @test_throws ArgumentError BandLimitedOneOverFNoiseModel(1.0, 0.1, 1.0, 0)
end

@testset "closed-form theory" begin
    Q = 0.7
    m = mathematica_model(Q)
    @test stationary_variance(m) ≈ Q * log(100)
    @test theoretical_autocovariance(m, 0.0) ≈ Q * log(100)
    @test theoretical_autocovariance(m, 1e-8) ≈ theoretical_autocovariance(m, 0.0) rtol = 1e-9  # continuous at 0
    for τ in (0.05, 0.3, 2.0)
        @test theoretical_autocovariance(m, τ) ≈ Q * (QNSb.cosint(100τ) - QNSb.cosint(τ))   # Mathematica kernel
        @test theoretical_autocovariance(m, -τ) == theoretical_autocovariance(m, τ)
        # Ci formula == ∫_{fl}^{fh} (Q/f) cos(2πfτ) df, evaluated in u = ln f
        u = range(log(m.fl), log(m.fh); length=200_001)
        @test trapz(u, [Q * cos(2π * exp(x) * τ) for x in u]) ≈ theoretical_autocovariance(m, τ) rtol = 1e-7
    end
    @test theoretical_psd(m, 1.0) ≈ Q / 1.0
    @test theoretical_psd(m, 0.1) == 0.0                                   # below fl
    @test theoretical_psd(m, 20.0) == 0.0                                  # above fh
    # ∫ S df = variance. Midpoint rule in u = ln f: never evaluates exactly at fl or fh,
    # where exp(log(fh)) can round to just above fh and the PSD correctly returns 0.
    n, h = 10_000, log(m.fh / m.fl) / 10_000
    @test sum(theoretical_psd(m, m.fl * exp((i - 0.5) * h)) * m.fl * exp((i - 0.5) * h) for i in 1:n) * h ≈ stationary_variance(m) rtol = 1e-8
end

@testset "synthesis grid" begin
    m = mathematica_model()
    n_t, dt = 1001, 1e-3
    g = QNSb._synthesis_grid(m, n_t, dt)
    df = 1 / (g.L * dt)
    @test g.L >= 8 * n_t
    @test sum(g.w) ≈ stationary_variance(m) rtol = 1e-12          # exact bin integrals ⇒ exact variance
    @test all(m.fl - df / 2 .<= g.f .<= m.fh + df / 2)
    @test all(g.w .> 0)
    # Quadrature error of the synthesised ACF vs the Ci formula, over lags up to the record length
    lags = (0:10:n_t-1) .* dt
    err = maximum(abs(generated_acov(g, τ) - theoretical_autocovariance(m, τ)) for τ in lags)
    @test err <= 5e-3 * stationary_variance(m)
    # Guards
    @test_throws ArgumentError QNSb._synthesis_grid(m, 101, 0.1)      # fh = 15.9 above Nyquist 5
    @test_throws ArgumentError QNSb._synthesis_grid(mathematica_model(; ext=1), 11, 1e-3)  # df ≫ fl
end

@testset "generated samples" begin
    m = mathematica_model(0.5)
    dt, M = 1e-3, 4000
    e = generate_ensemble(m, (0.0, 1.0), dt, M; rng=StableRNG(21))
    X = samples(e)
    @test size(X) == (1001, M)
    @test samples(generate_ensemble(m, (0.0, 1.0), dt, 3; rng=StableRNG(4))) ==
          samples(generate_ensemble(m, (0.0, 1.0), dt, 3; rng=StableRNG(4)))
    v = stationary_variance(m)
    for k in (1, 501, 1001)                                          # stationary everywhere
        @test within_se(mean(X[k, :]), 0.0, sqrt(v / M))
        @test within_se(var(X[k, :]), v, v * sqrt(2 / (M - 1)))
    end
    g = QNSb._synthesis_grid(m, 1001, dt)
    for lag in (30, 300)
        prods = X[101, :] .* X[101+lag, :]
        @test within_se(mean(prods), generated_acov(g, lag * dt), std(prods) / sqrt(M))
    end
end