# ─────────────────────────────────────────────────────────────────────────────
# What the validator reports: ⟨ρ_j(t_g)⟩, eigenvalues, and fidelity, with errors
# ─────────────────────────────────────────────────────────────────────────────
#
# These are the quantities of thesis §3.2.3–3.2.4, so they can be compared number
# for number with the teammates' master-equation results:
#   ⟨ρ_j(t_g)⟩ for the six cardinal states ............... Eq. 3.21
#   eigenvalues λ₁, λ₂ of ⟨ρ_j(t_g)⟩ .................... Figs. 3.2–3.19
#   average fidelity F̄ and error ⟨ε⟩_ξ = 1 − ⟨F̄⟩_ξ ...... Eqs. 3.24, 3.25, Figs. 3.20–3.22
# Each comes with a Monte-Carlo standard error. A master-equation result "agrees"
# only if it lies within a few SE; otherwise the difference is real.

"Apply a 4×4 superoperator to a 2×2 density matrix."
apply_channel(S::AbstractMatrix, ρ::AbstractMatrix) = unvec_dm(S * vec_dm(ρ))

"⟨S(t_g)⟩, the noise-averaged gate."
final_channel(r::GateResult) = r.S[end]

"""
    evolve_state(r::GateResult, ρ0) -> (times, ρ)

⟨ρ(t)⟩ at every saved time, in the frame of Eq. 3.12.
"""
evolve_state(r::GateResult, ρ0::AbstractMatrix) = (times=r.times, ρ=[apply_channel(S, ρ0) for S in r.S])

"""
    final_state(r::GateResult, ρ0) -> (ρ, se)

⟨ρ(t_g)⟩ = ⟨S(t_g)⟩[ρ0] (Eq. 3.21), together with the standard error of every
matrix element. `se[i, j]` is SE(Re ρᵢⱼ) + i·SE(Im ρᵢⱼ).
"""
function final_state(r::GateResult, ρ0::AbstractMatrix)
    v = vec_dm(ρ0)
    xs = [S * v for S in r.S_final]
    M = length(xs)
    m = sum(xs) / M
    se_re = sqrt.(sum(x -> abs2.(real.(x) .- real.(m)), xs) ./ (M - 1) ./ M)
    se_im = sqrt.(sum(x -> abs2.(imag.(x) .- imag.(m)), xs) ./ (M - 1) ./ M)
    return (ρ=unvec_dm(m), se=unvec_dm(complex.(se_re, se_im)))
end

"""
    fidelity_map(S, U_ideal) -> F̄                          (Eq. 3.24)

F̄ = (1/6) Σ_j Tr[U_ideal ρ_j U_ideal† ρ_j(t_g)], with ρ_j(t_g) = S[ρ_j], summed over
the six cardinal states. It works for ANY 4×4 map S: the brute-force average, one
trajectory, or a master-equation result.
"""
fidelity_map(S::AbstractMatrix, Uideal::AbstractMatrix) =
    sum(real(tr(Uideal * ρ * Uideal' * apply_channel(S, ρ))) for ρ in cardinal_states()) / 6

"""
    fidelity_trace_formula(S, U_ideal) = (Re Tr(S_ideal† S) + 2)/6

The same F̄, computed a second way through the entanglement fidelity (Horodecki et al.
1999; Nielsen 2002). The tests require it to equal [`fidelity_map`](@ref).
"""
fidelity_trace_formula(S::AbstractMatrix, Uideal::AbstractMatrix) =
    (real(tr(superoperator(Uideal)' * S)) + 2) / 6

"Every trajectory's own error ε_k = 1 − F̄(S_k(t_g)). Their mean is ⟨ε⟩_ξ."
trajectory_errors(r::GateResult) = [1 - fidelity_map(S, ideal_gate(r.gate)) for S in r.S_final]

"""
    average_error(r) -> (ε, se)        ⟨ε⟩_ξ = 1 − ⟨F̄⟩_ξ   (Eq. 3.25)
    average_fidelity(r) -> (F, se)     ⟨F̄⟩_ξ               (Eq. 3.24)

F̄ is linear in S, so averaging over trajectories commutes with Eq. 3.24. The mean
of the per-trajectory values is therefore exactly F̄ of the averaged channel, and
their spread gives the standard error.
"""
function average_error(r::GateResult)
    e = trajectory_errors(r)
    return (ε=mean(e), se=std(e) / sqrt(length(e)))
end

function average_fidelity(r::GateResult)
    e = average_error(r)
    return (F=1 - e.ε, se=e.se)
end

"""
    batch_estimate(f, r; nbatches=20) -> (value, se)

Standard error for a NONLINEAR function `f(S)` of the averaged channel, such as an
eigenvalue of ⟨ρ(t_g)⟩. The trajectories are split into `nbatches` groups, `f` is
evaluated on each group's average, and the spread of those values gives the SE
(the batch-means method). `value = f(full average)`. `f` must return a real number
or a real vector.
"""
function batch_estimate(f, r::GateResult; nbatches::Integer=20)
    M = length(r.S_final)
    2 <= nbatches <= M || throw(ArgumentError("need 2 ≤ nbatches ≤ number of trajectories ($M)"))
    edges = round.(Int, range(0, M; length=nbatches + 1))
    tovec(x) = x isa Number ? [Float64(x)] : Vector{Float64}(x)
    A = reduce(hcat, [tovec(f(sum(r.S_final[edges[b]+1:edges[b+1]]) / (edges[b+1] - edges[b])))
                      for b in 1:nbatches])
    se = vec(std(A; dims=2)) ./ sqrt(nbatches)
    v = f(final_channel(r))
    return v isa Number ? (value=v, se=se[1]) : (value=v, se=se)
end

"""
    final_state_eigenvalues(r, ρ0; nbatches=20) -> (λ, se)

Eigenvalues [λ₁, λ₂] of ⟨ρ(t_g)⟩ (as plotted in Figs. 3.2–3.19) and the SE of their real
parts. For the exact dynamics Im λ = 0 to rounding, and 0 ≤ λ ≤ 1.
"""
function final_state_eigenvalues(r::GateResult, ρ0::AbstractMatrix; nbatches::Integer=20)
    b = batch_estimate(S -> real.(density_matrix_eigenvalues(apply_channel(S, ρ0))), r; nbatches)
    return (λ=density_matrix_eigenvalues(apply_channel(final_channel(r), ρ0)), se=b.se)
end

"""
    choi_matrix(S) -> 4×4

J = Σᵢⱼ |i⟩⟨j| ⊗ S(|i⟩⟨j|). The map S is completely positive if and only if J ⪰ 0.
"""
function choi_matrix(S::AbstractMatrix)
    J = zeros(ComplexF64, 4, 4)
    for i in 1:2, j in 1:2
        E = zeros(ComplexF64, 2, 2)
        E[i, j] = 1
        J .+= kron(E, apply_channel(S, E))
    end
    return J
end

"""
    is_cptp(S; tol=1e-10) -> Bool

The map is completely positive (Choi matrix ⪰ 0) and trace preserving. The exact
average of unitary evolutions always passes. An approximate master equation may not
(thesis §3.3), and this is how to check it.
"""
function is_cptp(S::AbstractMatrix; tol::Real=1e-10)
    J = choi_matrix(S)
    maximum(abs, J - J') <= tol || return false
    cp = minimum(eigvals(Hermitian((J + J') / 2))) >= -tol
    tp = all(abs(tr(apply_channel(S, E)) - tr(E)) <= tol
             for E in (ComplexF64[1 0; 0 0], ComplexF64[0 1; 0 0], ComplexF64[0 0; 1 0], ComplexF64[0 0; 0 1]))
    return cp && tp
end