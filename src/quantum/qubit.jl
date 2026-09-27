# ─────────────────────────────────────────────────────────────────────────────
# Qubit basics: Pauli matrices, the six cardinal states, vectorisation, observables
# ─────────────────────────────────────────────────────────────────────────────

const σx = ComplexF64[0 1; 1 0]
const σy = ComplexF64[0 -im; im 0]
const σz = ComplexF64[1 0; 0 -1]
const I2 = ComplexF64[1 0; 0 1]

# ── Vectorisation (column-major, Julia's native order) ───────────────────────

"""
    vec_dm(ρ) -> Vector{ComplexF64}

Stack the columns of a 2×2 matrix: [ρ₁₁, ρ₂₁, ρ₁₂, ρ₂₂]. In this ordering
vec(AρB) = (Bᵀ ⊗ A)·vec(ρ). So ρ → UρU† is the 4×4 matrix conj(U) ⊗ U
([`superoperator`](@ref)), and the commutator ρ → [H, ρ] is 1 ⊗ H − Hᵀ ⊗ 1.
"""
vec_dm(ρ::AbstractMatrix) = vec(Matrix{ComplexF64}(ρ))

"Inverse of [`vec_dm`](@ref)."
function unvec_dm(v::AbstractVector)
    length(v) == 4 || throw(DimensionMismatch("expected a length-4 vector, got $(length(v))"))
    return reshape(Vector{ComplexF64}(v), 2, 2)
end

"`superoperator(U) = conj(U) ⊗ U`: the 4×4 matrix of ρ → UρU† acting on `vec_dm(ρ)`."
superoperator(U::AbstractMatrix) = kron(conj(U), U)

# ── States ───────────────────────────────────────────────────────────────────

"""
    cardinal_states() -> (zp, zm, xp, xm, yp, ym)

The six initial states of thesis §3.2.3 (Fig. 3.1), as density matrices, with the
thesis's names:

    ρ_zp = |0⟩⟨0|,                    ρ_zm = |1⟩⟨1|,
    ρ_xp = ½(|0⟩+|1⟩)(⟨0|+⟨1|),       ρ_xm = ½(|0⟩−|1⟩)(⟨0|−⟨1|),
    ρ_yp = ½(|0⟩+i|1⟩)(⟨0|−i⟨1|),     ρ_ym = ½(|0⟩−i|1⟩)(⟨0|+i⟨1|).

They form a 2-design: averaging a state fidelity over these six states equals
averaging it over the whole Bloch sphere. That is why Eq. 3.24 uses them.
"""
function cardinal_states()
    pure(a, b) = (ψ = ComplexF64[a, b]; ψ * ψ' / real(ψ' * ψ))
    return (zp=pure(1, 0), zm=pure(0, 1), xp=pure(1, 1), xm=pure(1, -1),
            yp=pure(1, im), ym=pure(1, -im))
end

"""
    bloch_state(θ, ϕ) -> ρ

Pure state cos(θ/2)|0⟩ + e^{iϕ} sin(θ/2)|1⟩ as a density matrix.
"""
function bloch_state(θ::Real, ϕ::Real)
    ψ = ComplexF64[cos(θ / 2), cis(ϕ) * sin(θ / 2)]
    return ψ * ψ'
end

# ── Observables ──────────────────────────────────────────────────────────────

"Population of |0⟩: ρ₁₁."
population_0(ρ) = real(ρ[1, 1])
"Population of |1⟩: ρ₂₂."
population_1(ρ) = real(ρ[2, 2])
"Purity Tr(ρ²): 1 for a pure state, ½ for the maximally mixed state."
purity(ρ) = real(tr(ρ * ρ))
"Bloch vector (⟨σx⟩, ⟨σy⟩, ⟨σz⟩)."
bloch_vector(ρ) = (real(tr(σx * ρ)), real(tr(σy * ρ)), real(tr(σz * ρ)))

"""
    density_matrix_eigenvalues(ρ) -> [λ₁, λ₂]   (complex, Re λ₁ ≥ Re λ₂)

Eigenvalues of ρ as a general matrix. It is deliberately NOT made Hermitian first,
so a density matrix from an approximate master equation shows the imaginary parts
the thesis plots (Figs. 3.2–3.19: Re λ₁, Re λ₂, Im λ₁, Im λ₂). A physical ρ has
λ₁, λ₂ ∈ [0, 1] with zero imaginary part.
"""
density_matrix_eigenvalues(ρ::AbstractMatrix) =
    sort(eigvals(Matrix{ComplexF64}(ρ)); by=real, rev=true)

"""
    validate_density_matrix(ρ; tol=1e-10) -> (hermitian_error, min_eigenvalue, trace_error, valid)

Checks the three physical requirements: ρ = ρ†, ρ ⪰ 0, and Tr ρ = 1.
"""
function validate_density_matrix(ρ::AbstractMatrix; tol::Real=1e-10)
    herm = maximum(abs, ρ - ρ')
    λmin = minimum(eigvals(Hermitian((ρ + ρ') / 2)))
    trerr = abs(tr(ρ) - 1)
    return (hermitian_error=herm, min_eigenvalue=λmin, trace_error=trerr,
            valid=herm <= tol && λmin >= -tol && trerr <= tol)
end