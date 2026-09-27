@testset "Pauli algebra" begin
    @test σx * σx == I2 && σy * σy == I2 && σz * σz == I2
    @test σx * σy ≈ im * σz
    @test σy * σz ≈ im * σx
    @test σz * σx ≈ im * σy
end

@testset "cardinal states are exactly the thesis §3.2.3 states" begin
    s = cardinal_states()
    @test keys(s) == (:zp, :zm, :xp, :xm, :yp, :ym)
    k0, k1 = ComplexF64[1, 0], ComplexF64[0, 1]
    @test s.zp ≈ k0 * k0'
    @test s.zm ≈ k1 * k1'
    @test s.xp ≈ (k0 + k1) * (k0 + k1)' / 2
    @test s.xm ≈ (k0 - k1) * (k0 - k1)' / 2
    @test s.yp ≈ (k0 + im * k1) * (k0 + im * k1)' / 2         # ½(|0⟩+i|1⟩)(⟨0|−i⟨1|)
    @test s.ym ≈ (k0 - im * k1) * (k0 - im * k1)' / 2
    blochs = ((0, 0, 1), (0, 0, -1), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))
    for (ρ, b) in zip(s, blochs)
        @test validate_density_matrix(ρ).valid
        @test purity(ρ) ≈ 1
        @test collect(bloch_vector(ρ)) ≈ collect(b) atol = 1e-15
    end
    @test bloch_state(π / 2, π / 2) ≈ s.yp
    @test population_0(s.zp) == 1 && population_1(s.zp) == 0
end

@testset "vectorisation is column-major" begin
    ρ = ComplexF64[1 2; 3 4]
    @test vec_dm(ρ) == [1, 3, 2, 4]
    @test unvec_dm(vec_dm(ρ)) == ρ
    U = su2_exp((0.3, -0.2, 0.7))
    ρ = bloch_state(1.1, 0.4)
    @test superoperator(U) * vec_dm(ρ) ≈ vec_dm(U * ρ * U')           # vec(UρU†) = (Ū ⊗ U) vec(ρ)
    H = ComplexF64[0.3 0.1-0.2im; 0.1+0.2im -0.5]
    @test (kron(I2, H) - kron(transpose(H), I2)) * vec_dm(ρ) ≈ vec_dm(H * ρ - ρ * H)   # commutator
    @test_throws DimensionMismatch unvec_dm(zeros(3))
end

@testset "eigenvalues keep imaginary parts; validity checks" begin
    @test density_matrix_eigenvalues(ComplexF64[0.8 0; 0 0.2]) ≈ [0.8, 0.2]
    λ = density_matrix_eigenvalues(ComplexF64[0.5 0.3; -0.3 0.5])    # not Hermitian
    @test imag(λ[1]) ≈ -imag(λ[2]) && abs(imag(λ[1])) ≈ 0.3
    @test !validate_density_matrix(ComplexF64[1.2 0; 0 -0.2]).valid   # negative eigenvalue
    @test !validate_density_matrix(ComplexF64[0.6 0; 0 0.6]).valid    # trace 1.2
    @test !validate_density_matrix(ComplexF64[0.5 0.1; 0.3 0.5]).valid # not Hermitian
end