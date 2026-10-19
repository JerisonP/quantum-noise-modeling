@testset "no noise: the ideal gate" begin
    g = CosineGate(π / 2, 1.0)
    me = tcl2_evolution(g, OUNoiseModel(0.0, 1.0); nsteps=200)
    @test me.times[end] ≈ 1.0 && length(me.S) == 201
    @test me.S[end] ≈ superoperator(ideal_gate(g))
end

@testset "θ = 0: the 2nd-order equation is exact for Gaussian noise" begin
    σ, τc = 2.0, 0.5
    me = tcl2_evolution(CosineGate(0.0, 1.0), OUNoiseModel(σ, τc); nsteps=1000)
    for (k, t) in ((501, 0.5), (1001, 1.0))
        exact = exp(-σ^2 * τc * (t - τc * (1 - exp(-t / τc)))) / 2
        @test me.times[k] ≈ t
        @test apply_channel(me.S[k], cardinal_states().xp)[1, 2] ≈ exact rtol = 1e-6
    end
end

@testset "numerically converged in nsteps; trace and Hermiticity preserved" begin
    g, m = CosineGate(π / 2, 1.0), OUNoiseModel(0.5, 1.0)
    ε(n) = 1 - fidelity_map(tcl2_evolution(g, m; nsteps=n).S[end], ideal_gate(g))
    @test abs(ε(500) / ε(2000) - 1) < 1e-4
    S = tcl2_evolution(g, m; nsteps=500).S[end]
    for ρ0 in cardinal_states()
        ρ = apply_channel(S, ρ0)
        @test tr(ρ) ≈ 1
        @test ρ ≈ ρ'
    end
end

@testset "weak noise: brute force and the 2nd-order equation agree" begin
    g, τc = CosineGate(π / 2, 1.0), 1.0
    δ = 0.02
    m = OUNoiseModel(sigma_for_delta(δ, τc, g.θ, g.tg), τc)
    r = simulate_gate(g, m, 2000; dt=1e-3, rng=StableRNG(31), save_every=1000)
    S2 = tcl2_evolution(g, m; nsteps=1000).S[end]
    mc = average_error(r)
    ε2 = 1 - fidelity_map(S2, ideal_gate(g))
    @test abs(mc.ε - ε2) <= 5mc.se + 0.05ε2                       # differences beyond 2nd order are O(δ²)
    for ρ0 in cardinal_states()
        fs = final_state(r, ρ0)
        @test all(abs.(real.(fs.ρ - apply_channel(S2, ρ0))) .<= 5 .* real.(fs.se) .+ 1e-4)
    end
end

@testset "TCL2 with nonzero-mean OU noise matches brute force" begin
    # A mean μ adds the first-order term −iμ[V(t), ρ] (a static detuning). Without it,
    # TCL2 is blind to μ and its error is −22% (μ = σ/2) and −54% (μ = σ) in ε.
    g, τc = CosineGate(π / 2, 1.0), 1.0
    σ = sigma_for_delta(0.02, τc, g.θ, g.tg)
    ε0 = 1 - fidelity_map(tcl2_evolution(g, OUNoiseModel(σ, τc); nsteps=1000).S[end], ideal_gate(g))
    for μ in (0.5σ, σ)
        m = OUNoiseModel(μ, σ, τc)
        r = simulate_gate(g, m, 4000; dt=2e-3, rng=StableRNG(41), save_every=500)
        S2 = tcl2_evolution(g, m; nsteps=1000).S[end]
        mc = average_error(r)
        ε2 = 1 - fidelity_map(S2, ideal_gate(g))
        @test abs(mc.ε - ε2) <= 5mc.se + 0.02ε2                    # 2%: the μ = 0 offset is ≈ 1.2%
        @test ε2 > 1.2ε0                                            # fails if TCL2 ignores μ (then ε2 == ε0)
        for ρ0 in cardinal_states()
            fs = final_state(r, ρ0)
            d = fs.ρ - apply_channel(S2, ρ0)
            @test all(abs.(real.(d)) .<= 5 .* real.(fs.se) .+ 1e-4)
            @test all(abs.(imag.(d)) .<= 5 .* imag.(fs.se) .+ 1e-4)  # the phase μ·t lives here
        end
    end
end