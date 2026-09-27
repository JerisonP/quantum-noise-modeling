@testset "δ, Eq. 3.20" begin
    @test small_noise_parameter(1.0, 1.0, π / 2, 1.0) ≈ 1 / sqrt(1 + π^2 / 4)
    @test small_noise_parameter(1.0, 1.0, π / 2, 1.0) ≈ 0.5370292721463151 rtol = 1e-14
    @test small_noise_parameter(0.3, 0.1, π / 2, 1.0) ≈ 0.008890980317281631 rtol = 1e-14
    @test small_noise_parameter(CosineGate(π / 2, 1.0), OUNoiseModel(0.3, 0.1)) ==
          small_noise_parameter(0.3, 0.1, π / 2, 1.0)
    for (δ, τc, tg) in ((0.1, 1.0, 1.0), (0.05, 1.0, 0.1), (0.003, 1.0, 0.01))
        @test small_noise_parameter(sigma_for_delta(δ, τc, π / 2, tg), τc, π / 2, tg) ≈ δ
    end
end

@testset "the integral behind δ (Eqs. 3.16–3.19), checked numerically" begin
    θ, tg, τc, t = π / 2, 1.0, 0.7, 0.9
    x = range(0, t; length=200_001)
    integral = trapz(x, [cis(θ * s / tg) * exp(-s / τc) for s in x])
    a = im * θ / tg - 1 / τc
    closed = (1 - exp(a * t)) * tg * τc / (tg - im * θ * τc)   # Eq. 3.17 (with its overall sign fixed)
    @test integral ≈ closed rtol = 1e-9
    D = tg^2 + θ^2 * τc^2                                       # Eq. 3.18 (with the cos-term sign fixed)
    @test abs2(closed) ≈ (-2exp(-t / τc) * tg^2 * τc^2 * cos(t * θ / tg) + 2exp(-t / τc) * tg^2 * τc^2 * cosh(t / τc)) / D
    t_long = 30τc                                               # Eq. 3.19: drop the terms decaying as e^{−t/τc}
    @test abs2((1 - exp(a * t_long)) * tg * τc / (tg - im * θ * τc)) ≈ tg^2 * τc^2 / D rtol = 1e-10
end

@testset "delta_sweep: reference curves" begin
    g = CosineGate(π / 2, 1.0)
    δs = [0.02, 0.1, 0.3]
    rows = delta_sweep(g, 1.0, δs, 300; dt=1e-2, rng=StableRNG(21))
    @test [r.δ for r in rows] == δs
    @test all(small_noise_parameter(r.σ, 1.0, π / 2, 1.0) ≈ r.δ for r in rows)
    @test issorted([r.ε for r in rows])                          # common random numbers: smooth in δ
    @test all(r.ε_se > 0 for r in rows)
    s = rows[2].states
    @test keys(s) == (:zp, :zm, :xp, :xm, :yp, :ym)
    @test all(validate_density_matrix(x.ρ).valid for x in s)
    @test all(length(x.λ) == 2 && length(x.λ_se) == 2 for x in s)
    # a row equals running the pieces by hand on the same scaled noise
    unit = generate_ensemble(OUNoiseModel(1.0, 1.0), (0.0, 1.0), 1e-2, 300; rng=StableRNG(21))
    r = simulate_gate(g, NoiseEnsemble(times(unit), rows[2].σ .* samples(unit)))
    @test average_error(r).ε ≈ rows[2].ε
    @test final_state(r, cardinal_states().yp).ρ ≈ s.yp.ρ
end