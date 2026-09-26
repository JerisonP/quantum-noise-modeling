@testset "Student t quantiles (reference: scipy.stats.t.ppf)" begin
    @test student_t_quantile(0.975, 10) ≈ 2.228138851986274 rtol = 1e-8
    @test student_t_quantile(0.975, 1) ≈ 12.706204736174694 rtol = 1e-8
    @test student_t_quantile(0.995, 5) ≈ 4.032142983555228 rtol = 1e-8
    @test student_t_quantile(0.025, 10) ≈ -2.2281388519862753 rtol = 1e-8
    @test student_t_quantile(0.975, 1000) ≈ 1.9623390808264083 rtol = 1e-8
    @test student_t_quantile(0.5, 3) == 0.0
    @test_throws ArgumentError student_t_quantile(1.0, 3)
end

@testset "ols (reference: scipy.stats.linregress)" begin
    x = [1.0, 2, 3, 4, 5, 6]
    y = [1.1, 2.3, 2.8, 4.5, 4.9, 6.4]
    r = ols(x, y)
    @test r.slope ≈ 1.0285714285714287 rtol = 1e-12
    @test r.intercept ≈ 0.06666666666666599 rtol = 1e-9
    @test r.se_slope ≈ 0.07358645246507342 rtol = 1e-10
    @test r.se_intercept ≈ 0.28657805939566167 rtol = 1e-10
    @test r.r2 ≈ 0.9799374936989618 rtol = 1e-12
    tq = student_t_quantile(0.975, 4)
    @test collect(r.slope_ci) ≈ [r.slope - tq * r.se_slope, r.slope + tq * r.se_slope]
    exact = ols(x, 2 .+ 3 .* x)
    @test all(isapprox.((exact.slope, exact.intercept, exact.r2), (3.0, 2.0, 1.0)))
    @test exact.se_slope < 1e-12
    @test_throws ArgumentError ols([1.0, 2.0], [1.0, 2.0])
end

@testset "powerlaw_fit" begin
    f = (1:2000) .* 0.01
    S = 3 .* f .^ -1.3
    raw = powerlaw_fit(f, S)
    @test raw.β ≈ 1.3 rtol = 1e-12
    @test 10^raw.log10_amplitude ≈ 3.0 rtol = 1e-10
    # Default band: central 40% of the log-frequency range
    lo, hi = log10(f[1]), log10(f[end])
    @test log10(raw.fmin) >= lo + 0.3 * (hi - lo) && log10(raw.fmax) <= hi - 0.3 * (hi - lo)
    # Log-binning averages S inside each bin before taking logs: a tiny,
    # known curvature bias (≈4e-5 here). This is why β̂ must be compared with the
    # same fit applied to the expected spectrum, not with α.
    @test powerlaw_fit(f, S; nbins=24).β ≈ 1.3 rtol = 1e-4
    band = powerlaw_fit(f, S; fmin=0.5, fmax=5.0)
    @test band.fmin >= 0.5 && band.fmax <= 5.0
    @test_throws ArgumentError powerlaw_fit(f, S; fmin=1.0, fmax=1.01)
end

@testset "exponential_fit" begin
    lags = collect(0:0.01:3)
    C = 2 .* exp.(-lags ./ 0.5)
    r = exponential_fit(lags, C)
    @test r.τ_c ≈ 0.5 rtol = 1e-10
    @test r.C0 ≈ 2.0 rtol = 1e-10
    @test r.n == count(<(0.5 * log(20)), lags)             # the leading run with C > 0.05·C(0)
    C2 = copy(C)
    C2[250] = 1.5                                          # a spurious tail value above threshold
    @test exponential_fit(lags, C2).n == r.n               # is not used
    @test_throws ArgumentError exponential_fit(lags, reverse(C))
end