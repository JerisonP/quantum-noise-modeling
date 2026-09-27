const tgrid = collect(range(0, 1; length=101))
const Id4 = Matrix{ComplexF64}(I, 4, 4)

"Linear interpolation of samples y on grid t at x (the same noise model as the solvers)."
function lin_interp(t, y, x)
    k = clamp(searchsortedlast(t, x), 1, length(t) - 1)
    return y[k] + (y[k+1] - y[k]) * (x - t[k]) / (t[k+1] - t[k])
end

"""
A third, deliberately simple solver: exponential midpoint rule, exp(−iH(t_mid)h), on a
fine grid, for ANY Hamiltonian function H(ξ, t). It exists so that the tests can swap in a
wrong Hamiltonian as a negative control.
"""
function midpoint_reference(H, t, ξ; n=20_000)
    tt = range(t[1], t[end]; length=n + 1)
    U = Matrix{ComplexF64}(I2)
    for k in 1:n
        tm = (tt[k] + tt[k+1]) / 2
        U = exp(-im * (tt[k+1] - tt[k]) * H(lin_interp(t, ξ, tm), tm)) * U
    end
    return superoperator(U)
end

@testset "su2_exp" begin
    v = (0.3, -1.2, 0.5)
    U = su2_exp(v)
    @test U ≈ exp(-im * (v[1] * σx + v[2] * σy + v[3] * σz))
    @test su2_exp((0.0, 0.0, 0.0)) == I2
end

@testset "zero noise: both solvers give exactly the ideal gate" begin
    g = CosineGate(π / 2, 1.0)
    for solver in (:magnus, :rk4)
        Ss = propagate_trajectory(g, tgrid, zeros(101); solver, substeps=4)
        @test Ss[1] ≈ Id4
        @test Ss[end] ≈ superoperator(ideal_gate(g)) atol = 1e-8
    end
end

@testset "no drive, constant detuning: exact precession" begin
    g0, ξ = CosineGate(0.0, 1.0), 2.7
    for solver in (:magnus, :rk4), (t, S) in zip(tgrid, propagate_trajectory(g0, tgrid, fill(ξ, 101); solver, substeps=8))
        @test S ≈ superoperator([cis(-ξ * t / 2) 0; 0 cis(ξ * t / 2)]) atol = 1e-9
    end
end

@testset "two independent solvers of Eq. 3.12 agree on every trajectory" begin
    g = CosineGate(π / 2, 1.0)
    rng = StableRNG(1)
    for _ in 1:3
        ξ = 3 .* randn(rng, 101)                                     # strong, rough noise
        Sm = propagate_trajectory(g, tgrid, ξ; solver=:magnus, substeps=8)[end]
        Sr = propagate_trajectory(g, tgrid, ξ; solver=:rk4, substeps=8)[end]
        @test Sm ≈ Sr atol = 1e-7
    end
end

@testset "the noise is on σz: agreement with a third solver, and a negative control" begin
    g = CosineGate(π / 2, 1.0)
    ξ = 3 .* randn(StableRNG(2), 101)
    S = propagate_trajectory(g, tgrid, ξ)[end]
    @test S ≈ midpoint_reference((z, t) -> hamiltonian(g, z, t), tgrid, ξ) atol = 2e-5
    on_σx = (z, t) -> drive(g, t) / 2 * σx + z / 2 * σx
    on_σy = (z, t) -> drive(g, t) / 2 * σx + z / 2 * σy
    @test maximum(abs, S - midpoint_reference(on_σx, tgrid, ξ)) > 0.05
    @test maximum(abs, S - midpoint_reference(on_σy, tgrid, ξ)) > 0.05
end

@testset "both solvers are 4th order" begin
    g = CosineGate(π / 2, 1.0)
    ξ = 3 .* randn(StableRNG(3), 101)
    ref = propagate_trajectory(g, tgrid, ξ; substeps=64)[end]
    for solver in (:magnus, :rk4)
        err(s) = maximum(abs, propagate_trajectory(g, tgrid, ξ; solver, substeps=s)[end] - ref)
        e1, e2, e4 = err(1), err(2), err(4)
        @test 3.7 < log2(e1 / e2) < 4.3 && 3.7 < log2(e2 / e4) < 4.3   # halving h cuts error 16×
    end
    @test_throws ArgumentError propagate_trajectory(g, tgrid, ξ; solver=:euler)
end

@testset "simulate_gate: grid checks, saved times, both solvers" begin
    g = CosineGate(π / 2, 1.0)
    ens = generate_ensemble(OUNoiseModel(1.0, 0.3), (0.0, 1.0), 1e-2, 50; rng=StableRNG(4))
    rm = simulate_gate(g, ens; save_every=10)
    rr = simulate_gate(g, ens; solver=:rk4, substeps=4, save_every=10)
    @test length(rm.times) == 11 && rm.times[end] ≈ 1.0
    @test rm.S[1] ≈ Id4
    @test length(rm.S_final) == 50 && rm.S_final[1] ≈ propagate_trajectory(g, times(ens), ens[1])[end]
    @test all(isapprox(a, b; atol=1e-6) for (a, b) in zip(rm.S, rr.S))
    @test final_channel(rm) ≈ sum(rm.S_final) / 50
    @test_throws ArgumentError simulate_gate(CosineGate(π / 2, 2.0), ens)     # grid ≠ [0, t_g]
end

@testset "exact OU dephasing (θ = 0): brute force reproduces the closed form" begin
    # No drive: ρ₀₁(T) = e^{−iΦ}ρ₀₁(0) with Φ = ∫₀ᵀ ξ. For Gaussian ξ,
    # ⟨ρ₀₁(T)⟩ = ½ exp(−Var Φ / 2), and for OU Var Φ = 2σ²τ_c (T − τ_c(1 − e^{−T/τ_c})).
    σ, τc, M = 2.0, 0.5, 4000
    for T in (0.5, 1.0)
        r = simulate_gate(CosineGate(0.0, T), OUNoiseModel(σ, τc), M; dt=0.01, rng=StableRNG(5))
        fs = final_state(r, cardinal_states().xp)
        exact = exp(-σ^2 * τc * (T - τc * (1 - exp(-T / τc)))) / 2
        @test within_se(real(fs.ρ[1, 2]), exact, real(fs.se[1, 2]))
        @test within_se(imag(fs.ρ[1, 2]), 0.0, imag(fs.se[1, 2]))  # ξ → −ξ symmetry
        @test fs.ρ ≈ apply_channel(final_channel(r), cardinal_states().xp)
    end
end

@testset "timestep_convergence: paired, and shrinking with dt" begin
    g = CosineGate(π / 2, 1.0)
    ens = generate_ensemble(OUNoiseModel(3.0, 0.02), (0.0, 1.0), 1 / 400, 1000; rng=StableRNG(6))
    rows = timestep_convergence(g, ens; factors=(1, 2, 4, 8))
    @test [r.factor for r in rows] == [1, 2, 4, 8]
    @test rows[1].Δε == 0 && rows[2].dt ≈ 2 / 400
    @test abs(rows[2].Δε) < abs(rows[4].Δε)                          # error grows with dt
    @test rows[4].Δε_se < rows[4].ε_se                               # pairing removes Monte-Carlo noise
    @test_throws ArgumentError timestep_convergence(g, ens; factors=(1, 3))
end