const gR = CosineGate(π / 2, 1.0)
const rR = simulate_gate(gR, OUNoiseModel(1.0, 0.3), 400; dt=1e-2, rng=StableRNG(11), save_every=10)

@testset "the brute-force average is a physical channel at every time" begin
    for S in rR.S
        @test is_cptp(S)
        @test apply_channel(S, I2 / 2) ≈ I2 / 2                        # unital: an average of unitaries
    end
    for ρ0 in cardinal_states(), ρ in evolve_state(rR, ρ0).ρ
        @test validate_density_matrix(ρ).valid
    end
end

@testset "is_cptp rejects unphysical maps" begin
    transpose_map = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]  # ρ → ρᵀ: positive but not CP
    @test !is_cptp(transpose_map)
    @test !is_cptp(2 .* Matrix{ComplexF64}(I, 4, 4))               # doubles the trace
    leaky = copy(final_channel(rR)); leaky[3, 1] += 0.01im          # breaks Hermiticity of ρ
    @test !is_cptp(leaky)
end

@testset "Eq. 3.24 fidelity: two formulas, per trajectory and averaged" begin
    Uid = ideal_gate(gR)
    S = final_channel(rR)
    @test fidelity_map(S, Uid) ≈ fidelity_trace_formula(S, Uid)      # 6 cardinal states = 2-design
    @test fidelity_map(superoperator(Uid), Uid) ≈ 1                  # perfect gate
    e = trajectory_errors(rR)
    @test all(x -> -1e-12 <= x <= 1, e)
    @test average_error(rR).ε ≈ 1 - fidelity_map(S, Uid)            # linear ⇒ mean commutes with Eq. 3.24
    @test average_error(rR).se ≈ std(e) / sqrt(400)
    @test average_fidelity(rR).F ≈ 1 - average_error(rR).ε
    # Eq. 3.24 written out in the thesis's own terms, state by state
    manual = sum(real(tr(Uid * ρ * Uid' * final_state(rR, ρ).ρ)) for ρ in cardinal_states()) / 6
    @test manual ≈ average_fidelity(rR).F
end

@testset "final states and eigenvalues (Eq. 3.21, Figs. 3.2–3.19)" begin
    for ρ0 in cardinal_states()
        fs = final_state(rR, ρ0)
        @test fs.ρ ≈ apply_channel(final_channel(rR), ρ0)
        @test all(real.(fs.se) .>= 0) && all(imag.(fs.se) .>= 0)
        ev = final_state_eigenvalues(rR, ρ0)
        @test sum(ev.λ) ≈ 1
        @test all(abs.(imag.(ev.λ)) .< 1e-12)                      # exact dynamics: no imaginary parts
        @test all(-1e-12 .<= real.(ev.λ) .<= 1 + 1e-12)
        @test real(ev.λ[1]) >= real(ev.λ[2])
        @test length(ev.se) == 2 && all(ev.se .> 0)
    end
end

@testset "batch_estimate" begin
    b = batch_estimate(S -> fidelity_map(S, ideal_gate(gR)), rR; nbatches=20)
    @test b.value ≈ average_fidelity(rR).F
    @test 0.3 < b.se / average_fidelity(rR).se < 3                  # batch SE ≈ direct SE for a linear quantity
    @test_throws ArgumentError batch_estimate(S -> 0.0, rR; nbatches=1)
end