const gP = CosineGate(π / 2, 1.0)
const mP = OUNoiseModel(1.0, 0.2)
const ensP = generate_ensemble(mP, (0.0, 1.0), 1e-2, 100; rng=StableRNG(61))

@testset "without Plots: a clear error, not a MethodError" begin
    @test Base.get_extension(QuantumNoiseSimulator, :QuantumNoiseSimulatorPlotsExt) === nothing
    err = try
        plot_noise_traces(ensP)
    catch e
        e
    end
    @test err isa ErrorException && occursin("using Plots", err.msg)
end

ENV["GKSwstype"] = "100"          # GR draws off-screen (CI machines have no display)
using Plots

@testset "loading Plots activates the extension" begin
    @test Base.get_extension(QuantumNoiseSimulator, :QuantumNoiseSimulatorPlotsExt) !== nothing
end

@testset "every plot builds, and one renders to a file" begin
    r = simulate_gate(gP, ensP; save_every=10)
    rows = delta_sweep(gP, 1.0, [0.02, 0.1, 0.2], 100; dt=1e-2, rng=StableRNG(62), nbatches=5)
    me = tcl2_delta_sweep(gP, 1.0, [0.02, 0.1, 0.2]; nsteps=200)
    plots = (plot_noise_traces(ensP; n_show=3),
             plot_validation(validate_noise_model(mP, ensP)),
             plot_state_evolution(r, cardinal_states().xp),
             plot_timestep_convergence(timestep_convergence(gP, ensP; factors=(1, 2, 4))),
             plot_error_vs_delta(rows; compare="2nd order" => me),
             plot_eigenvalues_vs_delta(rows, :xp; compare=["2nd order" => me]))
    for p in plots
        @test p isa Plots.Plot
    end
    @test length(plots[1].series_list) == 4                        # 3 traces + the zero line
    f = joinpath(mktempdir(), "error_vs_delta.png")
    savefig(plots[5], f)
    @test isfile(f) && filesize(f) > 1000
    @test_throws ArgumentError plot_eigenvalues_vs_delta(rows, :bogus)
end