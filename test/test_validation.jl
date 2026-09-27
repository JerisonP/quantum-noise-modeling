@testset "Bonferroni thresholds" begin
    @test QuantumNoiseSimulator._bonferroni_z(1, 1e-3) ≈ 3.2905267314919 rtol = 1e-8   # two-sided 99.9%
    @test QuantumNoiseSimulator._bonferroni_z(300, 1e-3) > QuantumNoiseSimulator._bonferroni_z(1, 1e-3)
end

@testset "correct generators pass" begin
    # α = 1e-5 per check keeps the chance of a false alarm across this whole
    # testset around 2e-4, while wrong models (next testset) still fail clearly.
    cases = ((OUNoiseModel(0.3, 1.0, 0.5), 0.01, 1001),
             (WhiteNoiseModel(0.7), 0.01, 501),
             (FractionalNoiseModel(1.0, 1e-3), 1e-3, 1001),
             (FractionalNoiseModel(0.5, 1e-3), 1e-3, 1001),
             (BandLimitedOneOverFNoiseModel(0.5, 1 / (2π), 100 / (2π)), 1e-3, 1001))
    for (i, (m, dt, N)) in enumerate(cases)
        e = generate_ensemble(m, (0.0, (N - 1) * dt), dt, 1000; rng=StableRNG(100 + i))
        r = validate_noise_model(m, e; α=1e-5)
        @test passed(r)
        @test length(r.checks) == 4
        @test length(r.acf.expected) == length(r.acf.C)
        @test length(r.psd.expected) == length(r.psd.S)
    end
end

@testset "the old per-trajectory demeaning also passes, once its bias is in the expectation" begin
    m = OUNoiseModel(1.0, 0.5)
    e = generate_ensemble(m, (0.0, 10.0), 0.01, 1000; rng=StableRNG(7))
    @test passed(validate_noise_model(m, e; demean=:trajectory, α=1e-5))
end

@testset "wrong models fail" begin
    ou = generate_ensemble(OUNoiseModel(1.0, 0.5), (0.0, 10.0), 0.01, 1000; rng=StableRNG(8))
    r = validate_noise_model(OUNoiseModel(1.0, 0.7), ou; α=1e-5)       # wrong τ_c
    @test !passed(r)
    @test !r.checks[3].passed                                            # the ACF check catches it

    pink = generate_ensemble(FractionalNoiseModel(1.0, 1e-3), (0.0, 1.0), 1e-3, 1000; rng=StableRNG(9))
    @test !passed(validate_noise_model(FractionalNoiseModel(1.0, 1.3e-3), pink; α=1e-5))   # wrong amplitude
end

@testset "report display" begin
    m = OUNoiseModel(1.0, 0.5)
    r = validate_noise_model(m, generate_ensemble(m, (0.0, 5.0), 0.01, 200; rng=StableRNG(3)))
    txt = sprint(show, MIME("text/plain"), r)
    @test occursin("Noise validation", txt) && occursin("overall:", txt)
    @test occursin("ValidationReport(OUNoiseModel", sprint(show, r))
end