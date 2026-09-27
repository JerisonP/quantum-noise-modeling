const ioDir = mktempdir()

@testset "noise ensembles: bit-for-bit round trip" begin
    ens = generate_ensemble(OUNoiseModel(1.0, 0.3), (0.0, 1.0), 0.01, 7; rng=StableRNG(51))
    f = save_ensemble(joinpath(ioDir, "sub", "noise.csv"), ens; meta=("model" => "OU(σ = 1, τ_c = 0.3)",))
    e2 = load_ensemble(f)
    @test times(e2) == times(ens) && samples(e2) == samples(ens)       # ==, not ≈
    t = read_table(f)
    @test t.meta["kind"] == "ensemble" && t.meta["n_trajectories"] == "7"
    @test t.meta["model"] == "OU(σ = 1, τ_c = 0.3)"
    @test t.cols[1:2] == ["t", "xi_1"] && size(t.data) == (101, 8)
    single = subensemble(ens, [3])
    @test samples(load_ensemble(save_ensemble(joinpath(ioDir, "one.csv"), single))) == samples(single)
end

@testset "channels: exact, including awkward floats" begin
    A = ComplexF64[0.1+1e-300im -0.0 Inf NaN; 1/3 nextfloat(0.0) -1e308 0; 1 2 3 4; 5 6 7 8im]
    f = save_channel(joinpath(ioDir, "odd.csv"), [0.0, 0.5], [A, 2A]; meta=("note" => "test",))
    c = load_channel(f)
    @test c.times == [0.0, 0.5]
    @test isequal(c.S[1], A) && isequal(c.S[2], 2A)                 # isequal: NaN == NaN, −0.0 ≠ 0.0
    @test c.meta["note"] == "test"

    g = CosineGate(π / 2, 1.0)
    r = simulate_gate(g, OUNoiseModel(1.0, 0.3), 20; dt=0.01, rng=StableRNG(52), save_every=25)
    c = load_channel(save_channel(joinpath(ioDir, "gate.csv"), r; all_times=true))
    @test c.times == r.times && c.S == r.S
    @test parse(Float64, c.meta["theta"]) == π / 2 && c.meta["solver"] == "magnus"
    @test parse(Int, c.meta["n_trajectories"]) == 20
    c1 = load_channel(save_channel(joinpath(ioDir, "gate_final.csv"), r))
    @test c1.times == [1.0] && c1.S == [final_channel(r)]
end

@testset "columns are read by name (the old load_propagator bug)" begin
    f = joinpath(ioDir, "gate.csv")
    t = read_table(f)
    perm = reverse(1:length(t.cols))                                # scramble the column order
    f2 = QuantumNoiseSimulator._write_table(joinpath(ioDir, "scrambled.csv"), "channel",
                                            (), t.cols[perm], t.data[:, perm])
    @test load_channel(f2).S == load_channel(f).S
end

@testset "δ sweeps: brute force, TCL2, and the hand-in format" begin
    g = CosineGate(π / 2, 1.0)
    rows = delta_sweep(g, 1.0, [0.02, 0.1], 60; dt=0.02, rng=StableRNG(53), nbatches=4)
    back = load_delta_sweep(save_delta_sweep(joinpath(ioDir, "sweep.csv"), rows; meta=("tau_c" => 1.0,)))
    @test length(back) == 2
    for (a, b) in zip(rows, back)
        @test (a.δ, a.σ, a.ε, a.ε_se) == (b.δ, b.σ, b.ε, b.ε_se)
        for s in keys(a.states)
            x, y = getproperty(a.states, s), getproperty(b.states, s)
            @test x.ρ == y.ρ && x.se == y.se && x.λ == y.λ && x.λ_se == y.λ_se
        end
    end
    me = tcl2_delta_sweep(g, 1.0, [0.02, 0.1]; nsteps=100)
    mb = load_delta_sweep(save_delta_sweep(joinpath(ioDir, "tcl2.csv"), me))
    @test [x.ε for x in mb] == [x.ε for x in me]
    @test mb[1].states.yp.ρ == me[1].states.yp.ρ
    @test isnan(mb[1].ε_se) && all(isnan, real.(mb[1].states.yp.se))  # a master equation has no SE

    # Hand-in format: only `delta` and `eps` are required
    f = joinpath(ioDir, "teammate.csv")
    write(f, """
    # QuantumNoiseSimulator CSV v1
    # kind = delta_sweep
    # source = SMNE 4th order
    delta,eps,xp_lambda1_re,xp_lambda2_re
    0.1,0.0165,0.97,0.03
    0.2,0.031,0.95,0.05
    """)
    tm = load_delta_sweep(f)
    @test [x.δ for x in tm] == [0.1, 0.2] && tm[2].ε == 0.031
    @test real(tm[1].states.xp.λ[1]) == 0.97 && isnan(imag(tm[1].states.xp.λ[1]))
    @test isnan(tm[1].σ) && all(isnan, real.(tm[1].states.zp.ρ))
end

@testset "bad files fail loudly" begin
    @test_throws ArgumentError load_channel(joinpath(ioDir, "sub", "noise.csv"))   # wrong kind
    plain = joinpath(ioDir, "plain.csv")
    write(plain, "a,b\n1,2\n")
    @test_throws ArgumentError read_table(plain)                                    # not our format
    noeps = joinpath(ioDir, "noeps.csv")
    write(noeps, "# QuantumNoiseSimulator CSV v1\n# kind = delta_sweep\ndelta\n0.1\n")
    @test_throws ArgumentError load_delta_sweep(noeps)                              # eps is required
    ragged = joinpath(ioDir, "ragged.csv")
    write(ragged, "# QuantumNoiseSimulator CSV v1\n# kind = delta_sweep\ndelta,eps\n0.1\n")
    @test_throws ArgumentError read_table(ragged)
end