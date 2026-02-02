"""
Superoperator Formalism for Density Matrix Evolution

This approach vectorizes the density matrix and computes an averaged Liouvillian
superoperator that can be applied to any initial state.
"""

using DifferentialEquations
using Plots
using LinearAlgebra
using Random
using Distributions
using Statistics
using Interpolations

# ============================================================================
# Vectorization and Matricization
# ============================================================================

"""
    vec_density_matrix(ρ::Matrix)

Vectorize a 2×2 density matrix into a 4×1 column vector.
Convention: |ρ⟩⟩ = [ρ00, ρ01, ρ10, ρ11]ᵀ
"""
function vec_density_matrix(ρ::Matrix{<:Number})
    return [ρ[1,1], ρ[1,2], ρ[2,1], ρ[2,2]]
end

"""
    unvec_density_matrix(ρ_vec::Vector)

Convert a vectorized density matrix back to 2×2 matrix form.
"""
function unvec_density_matrix(ρ_vec::Vector{<:Number})
    return [ρ_vec[1] ρ_vec[2]; ρ_vec[3] ρ_vec[4]]
end

"""
    vec_array_to_matrices(ρ_array)

Convert array of 8 real components to 2×2 complex matrix.
Input: [Re(ρ00), Im(ρ00), Re(ρ01), Im(ρ01), Re(ρ10), Im(ρ10), Re(ρ11), Im(ρ11)]
"""
function vec_array_to_matrix(u::Vector{<:Real})
    ρ = zeros(ComplexF64, 2, 2)
    ρ[1,1] = u[1] + im*u[2]  # ρ00
    ρ[1,2] = u[3] + im*u[4]  # ρ01
    ρ[2,1] = u[5] + im*u[6]  # ρ10
    ρ[2,2] = u[7] + im*u[8]  # ρ11
    return ρ
end

# ============================================================================
# Superoperator Construction
# ============================================================================

"""
    commutator_superoperator(H::Matrix)

Construct the superoperator for -i[H, ρ].

For a 2×2 Hamiltonian H, returns a 4×4 superoperator L_H such that:
vec(-i[H,ρ]) = L_H * vec(ρ)

Uses the identity: vec(AρB) = (Bᵀ ⊗ A)vec(ρ)
So: vec([H,ρ]) = vec(Hρ - ρH) = (I⊗H - Hᵀ⊗I)vec(ρ)
"""
function commutator_superoperator(H::Matrix{<:Number})
    I2 = Matrix{ComplexF64}(I, 2, 2)
    # vec([H,ρ]) = (I⊗H - Hᵀ⊗I)vec(ρ)
    # vec(-i[H,ρ]) = -i(I⊗H - Hᵀ⊗I)vec(ρ)
    L_comm = -im * (kron(I2, H) - kron(transpose(H), I2))
    return L_comm
end

"""
    build_liouvillian(H::Matrix)

Build the Liouvillian superoperator for coherent evolution.
L|ρ⟩⟩ = -i[H,ρ]|ρ⟩⟩
"""
function build_liouvillian(H::Matrix{<:Number})
    return commutator_superoperator(H)
end

"""
    build_time_dependent_liouvillian(t, tg, θ0)

Build the time-dependent Liouvillian for the driven qubit system.

H(t) = (ft/2) * σx where ft = (θ0/tg)*(1 - cos(2πt/tg))
"""
function build_time_dependent_liouvillian(t, tg, θ0)
    # Drive amplitude
    ft = (θ0/tg) * (1 - cos(2*π*t/tg))
    
    # Hamiltonian
    σx = [0.0+0im 1.0+0im; 1.0+0im 0.0+0im]
    H = (ft/2) * σx
    
    return build_liouvillian(H)
end

"""
    build_noisy_liouvillian(t, tg, θ0, ζ_t)

Build the Liouvillian with noise at time t.
Noise multiplicatively modifies the drive: ft → ft*(1 + ζ_t)
"""
function build_noisy_liouvillian(t, tg, θ0, ζ_t)
    # Drive amplitude with noise
    ft = (θ0/tg) * (1 - cos(2*π*t/tg))
    ft_noisy = ft * (1 + ζ_t)
    
    # Hamiltonian
    σx = [0.0+0im 1.0+0im; 1.0+0im 0.0+0im]
    H = (ft_noisy/2) * σx
    
    return build_liouvillian(H)
end

# ============================================================================
# Noise Generation
# ============================================================================

"""
    generate_ou_noise(τ_c, μ, σ, tspan, dt; kwargs...)

Generate Ornstein-Uhlenbeck noise trajectory.
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
    generate_ou_ensemble(τ_c, μ, σ, tspan, dt, n_samples; kwargs...)

Generate ensemble of OU noise trajectories.
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
# Average Liouvillian Computation
# ============================================================================

"""
    compute_average_liouvillian(times, tg, θ0, noise_ensemble)

Compute the time-dependent ensemble-averaged Liouvillian ⟨L(t)⟩.

# Arguments
- `times`: Time points at which to compute the average
- `tg`: Gate time
- `θ0`: Rotation angle
- `noise_ensemble`: EnsembleSolution of noise trajectories

# Returns
- Vector of 4×4 averaged Liouvillian matrices, one for each time point
"""
function compute_average_liouvillian(times, tg, θ0, noise_ensemble)
    n_times = length(times)
    n_samples = length(noise_ensemble)
    
    # Store average Liouvillian at each time
    L_avg = [zeros(ComplexF64, 4, 4) for _ in 1:n_times]
    
    println("Computing averaged Liouvillian over $n_samples noise realizations...")
    
    for i in 1:n_samples
        # Get noise trajectory
        noise_sol = noise_ensemble[i]
        noise_interp = linear_interpolation(noise_sol.t, noise_sol[1,:], 
                                           extrapolation_bc=Flat())
        
        # Compute Liouvillian for this noise realization at each time
        for (j, t) in enumerate(times)
            ζ_t = noise_interp(t)
            L_noisy = build_noisy_liouvillian(t, tg, θ0, ζ_t)
            L_avg[j] += L_noisy
        end
        
        if i % 10 == 0
            println("  Processed $i/$n_samples realizations")
        end
    end
    
    # Average
    for j in 1:n_times
        L_avg[j] /= n_samples
    end
    
    println("✓ Averaged Liouvillian computed!")
    
    return L_avg
end

# ============================================================================
# Evolution with Averaged Superoperator
# ============================================================================

"""
    evolve_with_averaged_superoperator(ρ0::Matrix, L_avg, times)

Evolve density matrix using the averaged superoperator.

# Arguments
- `ρ0`: Initial 2×2 density matrix
- `L_avg`: Vector of time-dependent averaged Liouvillians
- `times`: Time points

# Returns
- Array of evolved density matrices (one for each time point)
"""
function evolve_with_averaged_superoperator(ρ0::Matrix, L_avg, times)
    n_times = length(times)
    
    # Vectorize initial state (keep as ComplexF64)
    ρ_vec_0 = ComplexF64.(vec_density_matrix(ρ0))
    
    # Solve the ODE: dρ_vec/dt = L_avg(t) * ρ_vec
    function superop_ode!(dρ_vec, ρ_vec, p, t)
        times_ref, L_avg_ref = p
        # Find the Liouvillian at this time (simple nearest neighbor)
        t_idx = argmin(abs.(times_ref .- t))
        L_t = L_avg_ref[t_idx]
        
        # Compute derivative
        dρ_vec .= L_t * ρ_vec
        return nothing
    end
    
    # Solve
    tspan = (times[1], times[end])
    prob = ODEProblem(superop_ode!, ρ_vec_0, tspan, (times, L_avg))
    sol = solve(prob, Tsit5(), saveat=times)
    
    # Convert back to density matrices
    ρ_evolved = [unvec_density_matrix(sol.u[i]) for i in 1:length(sol.u)]
    
    return ρ_evolved, sol
end

# ============================================================================
# Standalone Superoperator Functions (What You Actually Want!)
# ============================================================================

"""
    compute_averaged_superoperator(tg, θ0, τ_c, σ, n_samples; n_times=1001)

Compute the time-dependent averaged Liouvillian superoperator ⟨L(t)⟩.

This is the MAIN function you want - it computes the superoperator once,
then you can apply it to any initial state!

# Arguments
- `tg`: Gate time
- `θ0`: Rotation angle  
- `τ_c`: Noise correlation time
- `σ`: Noise amplitude
- `n_samples`: Number of noise realizations to average over
- `n_times`: Number of time points

# Returns
Named tuple with:
- `L_avg`: Vector of 4×4 averaged Liouvillian matrices
- `times`: Time points
- `noise_ensemble`: The noise ensemble used (for reference)

# Example
```julia
# Compute superoperator once
sup = compute_averaged_superoperator(1.0, π/2, 0.1, 0.1, 100)

# Apply to different initial states
ρ_ground = apply_superoperator(ground_state_matrix(), sup)
ρ_excited = apply_superoperator(excited_state_matrix(), sup)
ρ_super = apply_superoperator(superposition_state_matrix(1.0, 1.0), sup)
```
"""
function compute_averaged_superoperator(tg, θ0, τ_c, σ, n_samples; n_times=1001)
    println("\n" * "="^60)
    println("Computing Averaged Liouvillian Superoperator")
    println("="^60)
    println("Parameters:")
    println("  tg = $tg, θ0 = $θ0")
    println("  τ_c = $τ_c, σ = $σ")  
    println("  n_samples = $n_samples")
    println("  n_times = $n_times")
    println()
    
    # Time grid
    times = collect(range(0, tg, length=n_times))
    dt = tg / 1000
    
    # Generate noise ensemble
    println("Step 1: Generating noise ensemble...")
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, (0.0, tg), dt, n_samples)
    println("✓ Noise ensemble generated!")
    println()
    
    # Compute averaged Liouvillian
    println("Step 2: Computing averaged Liouvillian...")
    L_avg = compute_average_liouvillian(times, tg, θ0, noise_ensemble)
    println()
    
    println("✓ Superoperator ready! Use apply_superoperator() to evolve states.")
    println()
    
    return (L_avg=L_avg, times=times, noise_ensemble=noise_ensemble)
end

"""
    apply_superoperator(ρ0::Matrix, superop)

Apply the averaged superoperator to evolve an initial density matrix.

# Arguments
- `ρ0`: Initial 2×2 density matrix
- `superop`: Result from compute_averaged_superoperator()

# Returns
Named tuple with:
- `ρ_evolved`: Array of density matrices at each time point
- `times`: Time points
- `sol`: Full ODE solution

# Example
```julia
sup = compute_averaged_superoperator(1.0, π/2, 0.1, 0.1, 50)
result = apply_superoperator(ground_state_matrix(), sup)
plot_superoperator_result(result)
```
"""
function apply_superoperator(ρ0::Matrix, superop)
    L_avg = superop.L_avg
    times = superop.times
    
    println("Evolving initial state using averaged superoperator...")
    ρ_evolved, sol = evolve_with_averaged_superoperator(ρ0, L_avg, times)
    println("✓ Evolution complete!")
    
    return (ρ_evolved=ρ_evolved, times=times, sol=sol, L_avg=L_avg)
end

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    superoperator_evolution(ρ0, tg, θ0, τ_c, σ, n_samples; n_times=1001)

Complete superoperator-based evolution with noise averaging.

# Steps:
1. Generate noise ensemble
2. Compute averaged Liouvillian ⟨L(t)⟩
3. Evolve any initial state ρ0 using ⟨L(t)⟩

# Arguments
- `ρ0`: Initial 2×2 density matrix
- `tg`: Gate time
- `θ0`: Rotation angle
- `τ_c`: Noise correlation time
- `σ`: Noise amplitude
- `n_samples`: Number of noise realizations
- `n_times`: Number of time points

# Returns
- Named tuple with (ρ_evolved, times, L_avg, sol)
"""
function superoperator_evolution(ρ0, tg, θ0, τ_c, σ, n_samples; n_times=1001)
    println("\n" * "="^60)
    println("Superoperator Formalism: Averaged Evolution")
    println("="^60)
    println("Parameters:")
    println("  tg = $tg, θ0 = $θ0")
    println("  τ_c = $τ_c, σ = $σ")
    println("  n_samples = $n_samples")
    println("  n_times = $n_times")
    println()
    
    # Time grid
    times = range(0, tg, length=n_times)
    dt = tg / 1000
    
    # Generate noise ensemble
    println("Step 1: Generating noise ensemble...")
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, (0.0, tg), dt, n_samples)
    println("✓ Noise ensemble generated!")
    println()
    
    # Compute averaged Liouvillian
    println("Step 2: Computing averaged Liouvillian...")
    L_avg = compute_average_liouvillian(times, tg, θ0, noise_ensemble)
    println()
    
    # Evolve with averaged superoperator
    println("Step 3: Evolving density matrix...")
    ρ_evolved, sol = evolve_with_averaged_superoperator(ρ0, L_avg, times)
    println("✓ Evolution complete!")
    println()
    
    return (ρ_evolved=ρ_evolved, times=collect(times), L_avg=L_avg, 
            sol=sol, noise_ensemble=noise_ensemble)
end

# ============================================================================
# Analysis and Plotting
# ============================================================================

"""
    extract_populations(ρ_evolved)

Extract populations from array of density matrices.
"""
function extract_populations(ρ_evolved)
    P0 = [real(ρ[1,1]) for ρ in ρ_evolved]
    P1 = [real(ρ[2,2]) for ρ in ρ_evolved]
    return (P0=P0, P1=P1)
end

"""
    plot_superoperator_result(result)

Plot the results from superoperator evolution.
"""
function plot_superoperator_result(result)
    pops = extract_populations(result.ρ_evolved)
    
    p = plot(result.times, pops.P0, label="ρ₀₀ (Ground)", 
             xlabel="Time", ylabel="Population",
             linewidth=2, title="Averaged Superoperator Evolution")
    plot!(p, result.times, pops.P1, label="ρ₁₁ (Excited)", linewidth=2)
    
    return p
end

"""
    compare_with_ensemble(ρ0, tg, θ0, τ_c, σ, n_samples; n_compare=10)

Compare superoperator method with direct ensemble simulation.
"""
function compare_with_ensemble(ρ0, tg, θ0, τ_c, σ, n_samples; n_compare=10)
    println("\n" * "="^60)
    println("Comparison: Superoperator vs Direct Ensemble")
    println("="^60)
    
    # Superoperator method
    println("\n--- Superoperator Method ---")
    result_super = superoperator_evolution(ρ0, tg, θ0, τ_c, σ, n_samples)
    pops_super = extract_populations(result_super.ρ_evolved)
    
    # Direct ensemble simulation (for comparison)
    println("\n--- Direct Ensemble Method ---")
    println("Simulating $n_compare trajectories for comparison...")
    
    # Load previous implementation functions
    include("density_matrix_simulation.jl")
    solutions_direct = simulate_noisy_ensemble(
        [real(ρ0[1,1]), imag(ρ0[1,1]), real(ρ0[1,2]), imag(ρ0[1,2]),
         real(ρ0[2,1]), imag(ρ0[2,1]), real(ρ0[2,2]), imag(ρ0[2,2])],
        tg, θ0, τ_c, σ, n_compare
    )
    stats_direct = ensemble_statistics(solutions_direct)
    
    # Plot comparison
    p = plot(result_super.times, pops_super.P0, 
             label="Superop: ρ₀₀", linewidth=3, linestyle=:dash,
             xlabel="Time", ylabel="Population",
             title="Superoperator vs Direct Ensemble")
    plot!(p, result_super.times, pops_super.P1, 
          label="Superop: ρ₁₁", linewidth=3, linestyle=:dash)
    
    plot!(p, stats_direct.times, stats_direct.P0_mean,
          ribbon=stats_direct.P0_std, fillalpha=0.2,
          label="Direct: ρ₀₀", linewidth=2, color=:blue)
    plot!(p, stats_direct.times, stats_direct.P1_mean,
          ribbon=stats_direct.P1_std, fillalpha=0.2,
          label="Direct: ρ₁₁", linewidth=2, color=:red)
    
    display(p)
    
    return (super=result_super, direct=stats_direct)
end

# ============================================================================
# Initial States
# ============================================================================

function ground_state_matrix()
    return ComplexF64[1.0 0.0; 0.0 0.0]
end

function excited_state_matrix()
    return ComplexF64[0.0 0.0; 0.0 1.0]
end

function superposition_state_matrix(α::Number, β::Number)
    norm = sqrt(abs2(α) + abs2(β))
    α, β = α/norm, β/norm
    ψ = [α; β]
    return ψ * ψ'  # |ψ⟩⟨ψ|
end

# ============================================================================
# Examples
# ============================================================================

"""
    example_superoperator()

Example: π/2 rotation with superoperator formalism.
"""
function example_superoperator()
    println("\n" * "="^60)
    println("Example: Superoperator Formalism (Recommended Workflow)")
    println("="^60)
    
    # Parameters
    tg = 1.0
    θ0 = π/2
    τ_c = 0.1
    σ = 0.1
    n_samples = 50
    
    # STEP 1: Compute superoperator ONCE
    println("\n*** STEP 1: Computing superoperator ***")
    sup = compute_averaged_superoperator(tg, θ0, τ_c, σ, n_samples, n_times=501)
    
    # STEP 2: Apply to initial state
    println("\n*** STEP 2: Applying to ground state ***")
    result = apply_superoperator(ground_state_matrix(), sup)
    
    # Final state
    ρ_final = result.ρ_evolved[end]
    println("\nFinal state:")
    println("  ρ₀₀ = $(ρ_final[1,1])")
    println("  ρ₁₁ = $(ρ_final[2,2])")
    println()
    
    # Plot
    display(plot_superoperator_result(result))
    
    # Return the superoperator so user can reuse it!
    return sup
end

"""
    example_reuse_superoperator()

Example showing how to reuse the same superoperator for multiple initial states.
This is the MAIN ADVANTAGE of the superoperator approach!
"""
function example_reuse_superoperator()
    println("\n" * "="^60)
    println("Example: Reusing Superoperator for Multiple States")
    println("="^60)
    
    # STEP 1: Compute superoperator ONCE  
    println("\n*** Computing superoperator ONCE ***")
    sup = compute_averaged_superoperator(1.0, π/2, 0.1, 0.1, 50, n_times=501)
    
    # STEP 2: Apply to MULTIPLE initial states
    println("\n*** Applying to 3 different initial states ***")
    
    println("\n  → Ground state...")
    result_ground = apply_superoperator(ground_state_matrix(), sup)
    
    println("\n  → Excited state...")
    result_excited = apply_superoperator(excited_state_matrix(), sup)
    
    println("\n  → Superposition state...")
    result_super = apply_superoperator(superposition_state_matrix(1.0, 1.0), sup)
    
    # Plot all three
    pops_ground = extract_populations(result_ground.ρ_evolved)
    pops_excited = extract_populations(result_excited.ρ_evolved)
    pops_super = extract_populations(result_super.ρ_evolved)
    
    p1 = plot(sup.times, pops_ground.P0, label="P₀", linewidth=2, 
              title="Initial: |0⟩", xlabel="Time", ylabel="Population")
    plot!(p1, sup.times, pops_ground.P1, label="P₁", linewidth=2)
    
    p2 = plot(sup.times, pops_excited.P0, label="P₀", linewidth=2,
              title="Initial: |1⟩", xlabel="Time", ylabel="Population")
    plot!(p2, sup.times, pops_excited.P1, label="P₁", linewidth=2)
    
    p3 = plot(sup.times, pops_super.P0, label="P₀", linewidth=2,
              title="Initial: (|0⟩+|1⟩)/√2", xlabel="Time", ylabel="Population")
    plot!(p3, sup.times, pops_super.P1, label="P₁", linewidth=2)
    
    display(plot(p1, p2, p3, layout=(3,1), size=(800, 900)))
    
    println("\n✓ All three states evolved using the SAME superoperator!")
    println("  This is MUCH faster than running 50 noise trajectories 3 times!")
    
    return sup
end

"""
    example_superoperator_old()

Example: π/2 rotation with superoperator formalism (old all-in-one version).
"""
function example_superoperator_old()
    println("\n" * "="^60)
    println("Example: Superoperator Formalism")
    println("="^60)
    
    # Parameters
    ρ0 = ground_state_matrix()
    tg = 1.0
    θ0 = π/2
    τ_c = 0.1
    σ = 0.1
    n_samples = 50
    
    # Run
    result = superoperator_evolution(ρ0, tg, θ0, τ_c, σ, n_samples, n_times=501)
    
    # Final state
    ρ_final = result.ρ_evolved[end]
    println("Final state:")
    println("  ρ₀₀ = $(ρ_final[1,1])")
    println("  ρ₁₁ = $(ρ_final[2,2])")
    println()
    
    # Plot
    display(plot_superoperator_result(result))
    
    return result
end

"""
    example_multiple_initial_states()

Demonstrate that the same averaged L can evolve different initial states.
"""
function example_multiple_initial_states()
    println("\n" * "="^60)
    println("Example: Multiple Initial States with Same ⟨L⟩")
    println("="^60)
    
    # Compute averaged Liouvillian once
    tg = 1.0
    θ0 = π/2
    τ_c = 0.1
    σ = 0.1
    n_samples = 50
    
    times = range(0, tg, length=501)
    dt = tg / 1000
    
    println("Computing averaged Liouvillian (once)...")
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, (0.0, tg), dt, n_samples)
    L_avg = compute_average_liouvillian(times, tg, θ0, noise_ensemble)
    
    # Try different initial states
    ρ0_ground = ground_state_matrix()
    ρ0_excited = excited_state_matrix()
    ρ0_super = superposition_state_matrix(1.0, 1.0)
    
    println("\nEvolving ground state...")
    ρ_ground, _ = evolve_with_averaged_superoperator(ρ0_ground, L_avg, collect(times))
    
    println("Evolving excited state...")
    ρ_excited, _ = evolve_with_averaged_superoperator(ρ0_excited, L_avg, collect(times))
    
    println("Evolving superposition state...")
    ρ_super, _ = evolve_with_averaged_superoperator(ρ0_super, L_avg, collect(times))
    
    # Extract and plot
    pops_ground = extract_populations(ρ_ground)
    pops_excited = extract_populations(ρ_excited)
    pops_super = extract_populations(ρ_super)
    
    p1 = plot(times, pops_ground.P0, label="P₀", linewidth=2, title="Initial: |0⟩")
    plot!(p1, times, pops_ground.P1, label="P₁", linewidth=2)
    
    p2 = plot(times, pops_excited.P0, label="P₀", linewidth=2, title="Initial: |1⟩")
    plot!(p2, times, pops_excited.P1, label="P₁", linewidth=2)
    
    p3 = plot(times, pops_super.P0, label="P₀", linewidth=2, title="Initial: (|0⟩+|1⟩)/√2")
    plot!(p3, times, pops_super.P1, label="P₁", linewidth=2)
    
    display(plot(p1, p2, p3, layout=(3,1), size=(800, 900)))
    
    println("\n✓ All three initial states evolved using the SAME averaged L!")
end

# ============================================================================
# Help Message
# ============================================================================

println("""
╔════════════════════════════════════════════════════════════╗
║  Superoperator Formalism Module Loaded!                   ║
╚════════════════════════════════════════════════════════════╝

RECOMMENDED WORKFLOW (Two-Step Process):

Step 1: Compute the averaged superoperator ONCE
    sup = compute_averaged_superoperator(tg, θ0, τ_c, σ, n_samples)
    
Step 2: Apply it to ANY initial states you want
    result1 = apply_superoperator(ground_state_matrix(), sup)
    result2 = apply_superoperator(excited_state_matrix(), sup)
    result3 = apply_superoperator(your_custom_state, sup)
    
This is much more efficient than running full ensemble for each state!

Plotting:
    plot_superoperator_result(result)
    
Alternative (All-in-one):
    result = superoperator_evolution(ρ0, tg, θ0, τ_c, σ, n_samples)
    
Examples:
    example_superoperator()              - Basic usage
    example_multiple_initial_states()    - Same ⟨L⟩ for different ρ0
    compare_with_ensemble()              - Compare methods

Key Functions:
    compute_averaged_superoperator()     - Get the superoperator
    apply_superoperator()                - Use it on any ρ0
""")