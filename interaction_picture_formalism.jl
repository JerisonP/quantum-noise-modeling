"""
Interaction Picture + Superoperator Formalism Library
WITH DRIVE AMPLITUDE NOISE (Multiplicative Noise Model)

This module provides tools for simulating noisy quantum gates using:
1. Interaction (rotating) picture transformation
2. Superoperator formalism with noise averaging
3. DRIVE AMPLITUDE NOISE: f(t) → f(t)*(1 + ζ(t))

Main workflow:
1. Compute averaged interaction picture superoperator: compute_interaction_picture_superoperator()
2. Apply to any initial state: apply_interaction_superoperator()

Author: Based on Mathematica code by Jerison L.
Date: 2025
"""

using DifferentialEquations
using LinearAlgebra
using Interpolations
using Statistics
using Distributions
using Random

# ============================================================================
# Constants and Basic Operators
# ============================================================================

const σx = ComplexF64[0 1; 1 0]
const σy = ComplexF64[0 -im; im 0]
const σz = ComplexF64[1 0; 0 -1]
const I2 = ComplexF64[1 0; 0 1]

# ============================================================================
# Interaction Picture Transformation Operators
# ============================================================================

"""
    U0(θ, tg, t)

Noiseless unitary evolution operator for interaction picture transformation.

For Hamiltonian H₀(t) = (1/2)(θ/tg)(1 - cos(2πt/tg))σx, the exact solution is:
    U₀(t) = cos(φ(t))I - i sin(φ(t))σx
where φ(t) = (θt)/(2tg) - (θ sin(2πt/tg))/(4π)

# Arguments
- `θ`: Rotation angle
- `tg`: Gate time
- `t`: Current time

# Returns
- 2×2 unitary matrix
"""
function U0(θ, tg, t)
    φ = (θ * t) / (2 * tg) - (θ * sin(2π * t / tg)) / (4π)
    return cos(φ) * I2 - im * sin(φ) * σx
end

"""
    U0_dagger(θ, tg, t)

Hermitian conjugate of U₀(t).
"""
function U0_dagger(θ, tg, t)
    φ = (θ * t) / (2 * tg) - (θ * sin(2π * t / tg)) / (4π)
    return cos(φ) * I2 + im * sin(φ) * σx
end

# ============================================================================
# Hamiltonians (WITH DRIVE AMPLITUDE NOISE)
# ============================================================================

"""
    H_lab(θ, tg, ζ, t)

Lab frame Hamiltonian with DRIVE AMPLITUDE NOISE.

H(t) = (1/2) * f(t) * (1 + ζ(t)) * σx

where f(t) = (θ/tg)(1 - cos(2πt/tg)) is the drive envelope
and ζ(t) is multiplicative noise on the drive amplitude.

# Arguments
- `θ`: Rotation angle
- `tg`: Gate time
- `ζ`: Multiplicative noise (ζ(t))
- `t`: Current time
"""
function H_lab(θ, tg, ζ, t)
    # Noiseless drive amplitude
    ft = (θ / tg) * (1 - cos(2π * t / tg))
    
    # Apply multiplicative noise
    ft_noisy = ft * (1 + ζ)
    
    return (ft_noisy / 2) * σx
end

"""
    H_interaction(θ, tg, ζ, t)

Interaction picture Hamiltonian with DRIVE AMPLITUDE NOISE.

H_I(t) = U₀†(t) [f(t)·ζ(t)/2 σx] U₀(t)

The noiseless drive is removed by U₀ transformation, leaving only
the noise-induced correction.

# Arguments
- `θ`: Rotation angle
- `tg`: Gate time
- `ζ`: Multiplicative noise ζ(t)
- `t`: Current time
"""
function H_interaction(θ, tg, ζ, t)
    U0_t = U0(θ, tg, t)
    U0d_t = U0_dagger(θ, tg, t)
    
    # Noiseless drive amplitude
    ft = (θ / tg) * (1 - cos(2π * t / tg))
    
    # Noise-induced term (the correction due to multiplicative noise)
    # H_noise = f(t) * ζ(t) / 2 * σx
    H_noise = (ft * ζ / 2) * σx
    
    # Transform to interaction picture
    return U0d_t * H_noise * U0_t
end

# ============================================================================
# Vectorization and Superoperator Tools
# ============================================================================

"""
    vec_density_matrix(ρ::Matrix)

Vectorize a 2×2 density matrix: |ρ⟩⟩ = [ρ₀₀, ρ₀₁, ρ₁₀, ρ₁₁]ᵀ
"""
function vec_density_matrix(ρ::Matrix{<:Number})
    return [ρ[1,1], ρ[1,2], ρ[2,1], ρ[2,2]]
end

"""
    unvec_density_matrix(ρ_vec::Vector)

Convert vectorized density matrix back to 2×2 matrix form.
"""
function unvec_density_matrix(ρ_vec::Vector{<:Number})
    return [ρ_vec[1] ρ_vec[2]; ρ_vec[3] ρ_vec[4]]
end

"""
    commutator_superoperator(H::Matrix)

Construct superoperator L_H such that vec(-i[H,ρ]) = L_H * vec(ρ).

Uses: vec([H,ρ]) = (I⊗H - Hᵀ⊗I)vec(ρ)
"""
function commutator_superoperator(H::Matrix{<:Number})
    I2_local = Matrix{ComplexF64}(I, 2, 2)
    L_comm = -im * (kron(I2_local, H) - kron(transpose(H), I2_local))
    return L_comm
end

"""
    build_liouvillian(H::Matrix)

Build Liouvillian superoperator: L|ρ⟩⟩ = -i[H,ρ]|ρ⟩⟩
"""
function build_liouvillian(H::Matrix{<:Number})
    return commutator_superoperator(H)
end

# ============================================================================
# Frame-Specific Liouvillians
# ============================================================================

"""
    build_interaction_liouvillian(t, tg, θ, ζ_t)

Build interaction picture Liouvillian at time t with noise ζ_t.

NOTE: Now using DRIVE AMPLITUDE NOISE (multiplicative)
"""
function build_interaction_liouvillian(t, tg, θ, ζ_t)
    H_I = H_interaction(θ, tg, ζ_t, t)
    return build_liouvillian(H_I)
end

"""
    build_lab_frame_liouvillian(t, tg, θ, ζ_t)

Build lab frame Liouvillian at time t with noise ζ_t.

NOTE: Now using DRIVE AMPLITUDE NOISE (multiplicative)
"""
function build_lab_frame_liouvillian(t, tg, θ, ζ_t)
    H_lab_t = H_lab(θ, tg, ζ_t, t)
    return build_liouvillian(H_lab_t)
end

# ============================================================================
# Noise Generation
# ============================================================================

"""
    generate_ou_noise(τ_c, μ, σ, tspan, dt; u0=0.0, solver=LambaEulerHeun(), 
                      reltol=1e-6, abstol=1e-8)

Generate single Ornstein-Uhlenbeck noise trajectory.

OU process: dX = (1/τc)(μ - X)dt + √(2σ²/τc)dW

This will be used as multiplicative noise on the drive amplitude.

# Arguments
- `τ_c`: Correlation time
- `μ`: Mean (typically 0)
- `σ`: Standard deviation
- `tspan`: Time span (t_start, t_end)
- `dt`: Sampling time step
"""
function generate_ou_noise(τ_c, μ, σ, tspan, dt; u0=0.0, solver=LambaEulerHeun(), 
                          reltol=1e-6, abstol=1e-8)
    f(u, p, t) = (1/τ_c) * (μ - u)
    g(u, p, t) = sqrt(2*σ*σ/τ_c)
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, saveat=dt, abstol=abstol, reltol=reltol)
    
    return sol
end

"""
    generate_ou_ensemble(τ_c, μ, σ, tspan, dt, n_samples; 
                         solver=LambaEulerHeun(), reltol=1e-6, abstol=1e-8)

Generate ensemble of OU noise trajectories with different initial conditions.

# Arguments
- `n_samples`: Number of independent noise realizations
"""
function generate_ou_ensemble(τ_c, μ, σ, tspan, dt, n_samples; 
                              solver=LambaEulerHeun(), reltol=1e-6, abstol=1e-8)
    f(u, p, t) = (1/τ_c) * (μ - u)
    g(u, p, t) = sqrt(2*σ*σ/τ_c)
    
    function prob_func(prob, i, repeat)
        remake(prob, u0=rand(Normal(0, σ)))
    end
    
    prob = SDEProblem(f, g, 0.0, tspan)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    
    sim = solve(ensemble_prob, solver, EnsembleThreads(), 
                trajectories=n_samples, saveat=dt, abstol=abstol, reltol=reltol)
    
    return sim
end

# ============================================================================
# Averaged Liouvillian Computation
# ============================================================================

"""
    compute_average_interaction_liouvillian(times, tg, θ, noise_ensemble)

Compute ensemble-averaged interaction picture Liouvillian ⟨L_I(t)⟩.

Uses DRIVE AMPLITUDE NOISE model: H = f(t)*(1 + ζ(t))/2 * σx

# Arguments
- `times`: Time points for evaluation
- `tg`: Gate time
- `θ`: Rotation angle
- `noise_ensemble`: EnsembleSolution of noise trajectories

# Returns
Vector of 4×4 averaged Liouvillian matrices at each time point
"""
function compute_average_interaction_liouvillian(times, tg, θ, noise_ensemble)
    n_times = length(times)
    n_samples = length(noise_ensemble)
    
    L_I_avg = [zeros(ComplexF64, 4, 4) for _ in 1:n_times]
    
    for i in 1:n_samples
        noise_sol = noise_ensemble[i]
        noise_interp = linear_interpolation(noise_sol.t, noise_sol[1,:], 
                                           extrapolation_bc=Flat())
        
        for (j, t) in enumerate(times)
            ζ_t = noise_interp(t)
            L_I_noisy = build_interaction_liouvillian(t, tg, θ, ζ_t)
            L_I_avg[j] += L_I_noisy
        end
    end
    
    for j in 1:n_times
        L_I_avg[j] /= n_samples
    end
    
    return L_I_avg
end

"""
    compute_average_lab_liouvillian(times, tg, θ, noise_ensemble)

Compute ensemble-averaged lab frame Liouvillian ⟨L_lab(t)⟩.

For comparison with interaction picture approach.
"""
function compute_average_lab_liouvillian(times, tg, θ, noise_ensemble)
    n_times = length(times)
    n_samples = length(noise_ensemble)
    
    L_lab_avg = [zeros(ComplexF64, 4, 4) for _ in 1:n_times]
    
    for i in 1:n_samples
        noise_sol = noise_ensemble[i]
        noise_interp = linear_interpolation(noise_sol.t, noise_sol[1,:], 
                                           extrapolation_bc=Flat())
        
        for (j, t) in enumerate(times)
            ζ_t = noise_interp(t)
            L_lab_noisy = build_lab_frame_liouvillian(t, tg, θ, ζ_t)
            L_lab_avg[j] += L_lab_noisy
        end
    end
    
    for j in 1:n_times
        L_lab_avg[j] /= n_samples
    end
    
    return L_lab_avg
end

# ============================================================================
# Evolution Functions
# ============================================================================

"""
    evolve_with_superoperator(ρ0_vec, L_avg, times)

Evolve vectorized density matrix using averaged superoperator.

Solves: dρ_vec/dt = L_avg(t) * ρ_vec

# Arguments
- `ρ0_vec`: Initial vectorized density matrix (4×1)
- `L_avg`: Vector of time-dependent Liouvillians
- `times`: Time points

# Returns
ODESolution object
"""
function evolve_with_superoperator(ρ0_vec, L_avg, times)
    function superop_ode!(dρ_vec, ρ_vec, p, t)
        times_ref, L_avg_ref = p
        t_idx = argmin(abs.(times_ref .- t))
        L_t = L_avg_ref[t_idx]
        dρ_vec .= L_t * ρ_vec
        return nothing
    end
    
    tspan = (times[1], times[end])
    prob = ODEProblem(superop_ode!, ρ0_vec, tspan, (times, L_avg))
    sol = solve(prob, Tsit5(), saveat=times)
    
    return sol
end

"""
    evolve_interaction_picture(ρ0_lab::Matrix, L_I_avg, times, tg, θ)

Evolve density matrix using interaction picture superoperator.

Steps:
1. Transform to interaction picture: ρ_I(0) = U₀†(0) ρ_lab(0) U₀(0) = ρ_lab(0)
2. Evolve: dρ_I/dt = ⟨L_I(t)⟩ ρ_I
3. Transform back: ρ_lab(t) = U₀(t) ρ_I(t) U₀†(t)

# Returns
Named tuple: (ρ_lab, ρ_I, times)
"""
function evolve_interaction_picture(ρ0_lab::Matrix, L_I_avg, times, tg, θ)
    # Step 1: Initial state (at t=0, U₀=I, so ρ_I(0) = ρ_lab(0))
    ρ_I_0_vec = ComplexF64.(vec_density_matrix(ρ0_lab))
    
    # Step 2: Evolve in interaction picture
    sol = evolve_with_superoperator(ρ_I_0_vec, L_I_avg, times)
    
    # Convert to matrices
    ρ_I_evolved = [unvec_density_matrix(sol.u[i]) for i in 1:length(sol.u)]
    
    # Step 3: Transform back to lab frame
    ρ_lab_evolved = similar(ρ_I_evolved)
    for (i, t) in enumerate(times)
        U0_t = U0(θ, tg, t)
        U0d_t = U0_dagger(θ, tg, t)
        ρ_lab_evolved[i] = U0_t * ρ_I_evolved[i] * U0d_t
    end
    
    return (ρ_lab=ρ_lab_evolved, ρ_I=ρ_I_evolved, times=collect(times))
end

"""
    evolve_lab_frame(ρ0_lab::Matrix, L_lab_avg, times)

Evolve density matrix using lab frame superoperator (for comparison).

# Returns
Named tuple: (ρ_lab, times)
"""
function evolve_lab_frame(ρ0_lab::Matrix, L_lab_avg, times)
    ρ0_vec = ComplexF64.(vec_density_matrix(ρ0_lab))
    sol = evolve_with_superoperator(ρ0_vec, L_lab_avg, times)
    ρ_lab_evolved = [unvec_density_matrix(sol.u[i]) for i in 1:length(sol.u)]
    return (ρ_lab=ρ_lab_evolved, times=collect(times))
end

# ============================================================================
# Main Interface Functions
# ============================================================================

"""
    compute_interaction_picture_superoperator(tg, θ, τ_c, σ, n_samples; 
                                              n_times=1001, verbose=true)

Compute averaged interaction picture superoperator.

This is the main function for the interaction picture approach.
Computes ⟨L_I(t)⟩ which can be reused for any initial state.

NOTE: Uses DRIVE AMPLITUDE NOISE model (multiplicative noise)

# Arguments
- `tg`: Gate time
- `θ`: Rotation angle  
- `τ_c`: Noise correlation time
- `σ`: Noise amplitude
- `n_samples`: Number of noise realizations
- `n_times`: Number of time points
- `verbose`: Print progress messages

# Returns
Named tuple with:
- `L_I_avg`: Averaged interaction picture Liouvillian
- `times`: Time grid
- `tg, θ`: System parameters (stored for back-transformation)
- `noise_ensemble`: Noise ensemble used

# Example
```julia
sup = compute_interaction_picture_superoperator(1.0, π/2, 0.1, 0.1, 100)
result = apply_interaction_superoperator(ground_state(), sup)
```
"""
function compute_interaction_picture_superoperator(tg, θ, τ_c, σ, n_samples; 
                                                   n_times=1001, verbose=true)
    verbose && println("="^60)
    verbose && println("Computing Interaction Picture Superoperator")
    verbose && println("="^60)
    verbose && println("Parameters: tg=$tg, θ=$θ")
    verbose && println("Noise: τ_c=$τ_c, σ=$σ (DRIVE AMPLITUDE)")
    verbose && println("Samples: $n_samples, Time points: $n_times")
    verbose && println()
    
    times = collect(range(0, tg, length=n_times))
    dt = tg / 1000
    
    verbose && println("Generating noise ensemble...")
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, (0.0, tg), dt, n_samples)
    
    verbose && println("Computing averaged interaction picture Liouvillian...")
    L_I_avg = compute_average_interaction_liouvillian(times, tg, θ, noise_ensemble)
    
    verbose && println("✓ Superoperator ready!")
    verbose && println()
    
    return (L_I_avg=L_I_avg, times=times, tg=tg, θ=θ, 
            noise_ensemble=noise_ensemble)
end

"""
    compute_lab_frame_superoperator(tg, θ, τ_c, σ, n_samples; 
                                    n_times=1001, verbose=true)

Compute averaged lab frame superoperator (for comparison).
"""
function compute_lab_frame_superoperator(tg, θ, τ_c, σ, n_samples; 
                                        n_times=1001, verbose=true)
    verbose && println("="^60)
    verbose && println("Computing Lab Frame Superoperator")
    verbose && println("="^60)
    
    times = collect(range(0, tg, length=n_times))
    dt = tg / 1000
    
    verbose && println("Generating noise ensemble...")
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, (0.0, tg), dt, n_samples)
    
    verbose && println("Computing averaged lab frame Liouvillian...")
    L_lab_avg = compute_average_lab_liouvillian(times, tg, θ, noise_ensemble)
    
    verbose && println("✓ Superoperator ready!")
    verbose && println()
    
    return (L_lab_avg=L_lab_avg, times=times, tg=tg, θ=θ,
            noise_ensemble=noise_ensemble)
end

"""
    apply_interaction_superoperator(ρ0_lab::Matrix, superop; verbose=true)

Apply interaction picture superoperator to evolve initial state.

# Arguments
- `ρ0_lab`: Initial 2×2 density matrix in lab frame
- `superop`: Output from compute_interaction_picture_superoperator()

# Returns
Named tuple: (ρ_lab, ρ_I, times)
"""
function apply_interaction_superoperator(ρ0_lab::Matrix, superop; verbose=true)
    verbose && println("Evolving state in interaction picture...")
    
    result = evolve_interaction_picture(ρ0_lab, superop.L_I_avg, 
                                       superop.times, superop.tg, superop.θ)
    
    verbose && println("✓ Evolution complete!")
    return result
end

"""
    apply_lab_superoperator(ρ0_lab::Matrix, superop; verbose=true)

Apply lab frame superoperator to evolve initial state.
"""
function apply_lab_superoperator(ρ0_lab::Matrix, superop; verbose=true)
    verbose && println("Evolving state in lab frame...")
    
    result = evolve_lab_frame(ρ0_lab, superop.L_lab_avg, superop.times)
    
    verbose && println("✓ Evolution complete!")
    return result
end

# ============================================================================
# Analysis Tools
# ============================================================================

"""
    extract_populations(ρ_evolved)

Extract ground and excited state populations from density matrices.

# Returns
Named tuple: (P0, P1) where P0 = ⟨0|ρ|0⟩ and P1 = ⟨1|ρ|1⟩
"""
function extract_populations(ρ_evolved)
    P0 = [real(ρ[1,1]) for ρ in ρ_evolved]
    P1 = [real(ρ[2,2]) for ρ in ρ_evolved]
    return (P0=P0, P1=P1)
end

"""
    extract_coherences(ρ_evolved)

Extract off-diagonal coherences from density matrices.

# Returns
Named tuple: (ρ01, ρ10)
"""
function extract_coherences(ρ_evolved)
    ρ01 = [ρ[1,2] for ρ in ρ_evolved]
    ρ10 = [ρ[2,1] for ρ in ρ_evolved]
    return (ρ01=ρ01, ρ10=ρ10)
end

"""
    purity(ρ::Matrix)

Compute purity Tr(ρ²) of a density matrix.
"""
function purity(ρ::Matrix)
    return real(tr(ρ * ρ))
end

"""
    fidelity(ρ1::Matrix, ρ2::Matrix)

Compute fidelity F(ρ1, ρ2) = Tr(√(√ρ1 ρ2 √ρ1)).

For simplicity, uses F = Tr(ρ1 ρ2) approximation for nearly pure states.
"""
function fidelity(ρ1::Matrix, ρ2::Matrix)
    return real(tr(ρ1 * ρ2))
end

# ============================================================================
# Initial States
# ============================================================================

"""
    ground_state()

Ground state density matrix |0⟩⟨0|
"""
ground_state() = ComplexF64[1.0 0.0; 0.0 0.0]

"""
    excited_state()

Excited state density matrix |1⟩⟨1|
"""
excited_state() = ComplexF64[0.0 0.0; 0.0 1.0]

"""
    superposition_state(α, β)

Superposition state density matrix |ψ⟩⟨ψ| where |ψ⟩ = α|0⟩ + β|1⟩
"""
function superposition_state(α::Number, β::Number)
    norm = sqrt(abs2(α) + abs2(β))
    α, β = α/norm, β/norm
    ψ = [α; β]
    return ψ * ψ'
end

"""
    thermal_state(T, ω)

Thermal state at temperature T with level spacing ω.
ρ_thermal = exp(-ω/T) / Z where Z = 1 + exp(-ω/T)
"""
function thermal_state(T, ω)
    if T == 0
        return ground_state()
    end
    Z = 1 + exp(-ω/T)
    p0 = 1/Z
    p1 = exp(-ω/T)/Z
    return ComplexF64[p0 0.0; 0.0 p1]
end

# ============================================================================
# Export Main Interface
# ============================================================================

export compute_interaction_picture_superoperator, compute_lab_frame_superoperator
export apply_interaction_superoperator, apply_lab_superoperator
export extract_populations, extract_coherences, purity, fidelity
export ground_state, excited_state, superposition_state, thermal_state
export U0, U0_dagger, H_lab, H_interaction
export generate_ou_noise, generate_ou_ensemble

println("✓ InteractionPictureSuperoperator.jl loaded successfully")