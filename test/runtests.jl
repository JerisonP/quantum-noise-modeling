using QuantumNoiseSimulator
using Test, Statistics, StableRNGs

# Statistical assertions in this suite follow one rule: the tolerance is k
# standard errors of the estimator (default k = 5), where the SE is computed
# from theory or from i.i.d. data. With a fixed StableRNG seed every test is
# deterministic; 5 SE also means an honest implementation would essentially
# never fail even if the seed changed (two-sided p ≈ 6e-7 per check).
"`est` agrees with `target` to within `k` standard errors `se`."
within_se(est, target, se; k=5) = abs(est - target) <= k * se

"Trapezoidal rule on a sorted grid."
trapz(x, y) = sum((x[i+1] - x[i]) * (y[i+1] + y[i]) / 2 for i in 1:length(x)-1)

@testset "QuantumNoiseSimulator" begin
    @testset "NoiseEnsemble" begin include("test_ensemble.jl") end
    @testset "OUNoiseModel" begin include("test_ou.jl") end
    @testset "WhiteNoiseModel" begin include("test_white.jl") end
    @testset "FractionalNoiseModel" begin include("test_fractional.jl") end
    @testset "BandLimitedOneOverFNoiseModel" begin include("test_bandlimited.jl") end
    @testset "Estimators" begin include("test_estimators.jl") end
    @testset "Fits" begin include("test_fits.jl") end
    @testset "Expectations" begin include("test_expectations.jl") end
    @testset "Validation" begin include("test_validation.jl") end
    @testset "Code quality (Aqua)" begin include("test_aqua.jl") end
end