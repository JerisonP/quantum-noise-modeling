const gH = CosineGate(π / 2, 1.0)            # the thesis gate: θ(t_g) = π/2

@testset "constructor" begin
    @test (gH.θ, gH.tg) == (π / 2, 1.0)
    @test_throws ArgumentError CosineGate(π, 0.0)
end

@testset "drive envelope and rotation angle (Eqs. 3.7, 3.15)" begin
    g = CosineGate(π / 2, 0.37)
    t = range(0, g.tg; length=100_001)
    @test trapz(t, [drive(g, x) for x in t]) ≈ g.θ rtol = 1e-8        # θ(t_g) = ∫ f_x = θ
    @test drive(g, 0.0) == 0 && abs(drive(g, g.tg)) < 1e-14            # turns on and off
    for x in (0.05, 0.2, 0.33)
        χ = g.θ / (2π) * sin(2π * x / g.tg)                           # Eq. 3.15
        @test rotation_angle(g, x) ≈ g.θ * x / g.tg - χ                # θ t/t_g − χ(t)
        h = 1e-7
        @test (rotation_angle(g, x + h) - rotation_angle(g, x - h)) / (2h) ≈ drive(g, x) rtol = 1e-7
    end
    @test rotation_angle(g, g.tg) ≈ g.θ
end

@testset "noiseless propagator (Eqs. 3.4–3.6, 3.8)" begin
    g = gH
    for x in (0.0, 0.25, 0.6, 1.0)
        U = ideal_propagator(g, x)
        @test U' * U ≈ I2
        @test U ≈ cos(rotation_angle(g, x) / 2) * I2 - im * sin(rotation_angle(g, x) / 2) * σx   # Eq. 3.6
        h = 1e-6
        dU = (ideal_propagator(g, x + h) - ideal_propagator(g, x - h)) / (2h)
        @test im * dU ≈ drive(g, x) / 2 * σx * U atol = 1e-8        # Eq. 3.5 with H_q,0 = f_x/2 σx (Eq. 3.4)
    end
    @test ideal_gate(g) ≈ (I2 - im * σx) / sqrt(2)                     # Eq. 3.8
end

@testset "the model, Eq. 3.12: the noise ξ sits on σz/2" begin
    g = gH
    ξ, t = 0.83, 0.41
    f = drive(g, t)
    @test hamiltonian(g, ξ, t) ≈ [ξ/2 f/2; f/2 -ξ/2]                  # f_x/2 σx + ξ/2 σz
    @test hamiltonian(g, 0.0, t) ≈ f / 2 * σx                          # H_q,0, Eq. 3.4
    @test hamiltonian(g, ξ, t) - hamiltonian(g, 0.0, t) ≈ ξ / 2 * σz   # V(t) = ½ξσz
end

@testset "interaction picture: U₀† V U₀, and the Eq. 3.13 ordering" begin
    g = gH
    ξ, t = 0.83, 0.41
    U0 = ideal_propagator(g, t)
    V = ξ / 2 * σz
    @test interaction_noise(g, ξ, t) ≈ U0' * V * U0                   # correct for i∂U₀ = H₀U₀
    a = rotation_angle(g, t)
    n = noise_axis(g, t)
    @test collect(n) ≈ [0, sin(a), cos(a)]
    @test U0' * σz * U0 ≈ n[2] * σy + n[3] * σz
    @test U0 * σz * U0' ≈ -n[2] * σy + n[3] * σz                       # Eq. 3.13 as printed: n_y flips sign
    @test noise_axis(CosineGate(0.0, 1.0), t) == (0.0, 0.0, 1.0)     # no drive: noise stays along z
end