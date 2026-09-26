@testset "construction and accessors" begin
    t = collect(0.0:0.1:1.0)
    X = reshape(collect(1.0:33.0), 11, 3)
    e = NoiseEnsemble(t, X)
    @test ntimes(e) == 11
    @test ntrajectories(e) == 3
    @test length(e) == 3
    @test timestep(e) ≈ 0.1
    @test times(e) == t
    @test samples(e) == X
    @test e[2] == X[:, 2]
    @test collect(e) == [X[:, j] for j in 1:3]          # iteration = columns
    @test sprint(show, e) == "NoiseEnsemble(3 trajectories × 11 time points, dt = 0.1)"
end

@testset "subensemble copies selected columns" begin
    e = NoiseEnsemble(0.0:1.0:4.0, reshape(collect(1.0:20.0), 5, 4))
    s = subensemble(e, [4, 1])
    @test samples(s) == samples(e)[:, [4, 1]]
    @test times(s) == times(e)
    samples(s)[1, 1] = -1.0
    @test samples(e)[1, 4] != -1.0                      # independent copy
end

@testset "invalid input is rejected" begin
    @test_throws DimensionMismatch NoiseEnsemble(0.0:1.0:3.0, zeros(5, 2))
    @test_throws ArgumentError NoiseEnsemble([0.0], zeros(1, 2))              # 1 point
    @test_throws ArgumentError NoiseEnsemble(0.0:1.0:3.0, zeros(4, 0))        # 0 trajectories
    @test_throws ArgumentError NoiseEnsemble([0.0, 0.1, 0.3, 0.4], zeros(4, 1))  # non-uniform
    @test_throws ArgumentError NoiseEnsemble([1.0, 0.0], zeros(2, 1))         # decreasing
end

@testset "time grid helper" begin
    tg = QuantumNoiseSimulator._time_grid((0.0, 1.0), 0.25)
    @test tg == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test length(QuantumNoiseSimulator._time_grid((0.0, 4.0), 1e-3)) == 4001
    @test_throws ArgumentError QuantumNoiseSimulator._time_grid((0.0, 1.0), 0.3)  # 3.33 steps
    @test_throws ArgumentError QuantumNoiseSimulator._time_grid((0.0, 1.0), -0.1)
    @test_throws ArgumentError QuantumNoiseSimulator._time_grid((1.0, 0.0), 0.1)
end