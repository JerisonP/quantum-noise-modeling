"""
Interaction Picture (Dirac Picture) Formalism

Splits Hamiltonian into H = H₀ + H₁(t):
- H₀: Free part (detuning) - transforms away
- H₁: Interaction part (drive) - treated perturbatively

Key transformation:
    ρ_I(t) = U₀†(t) ρ_S(t) U₀(t)  where U₀ = exp(-iH₀t)
    
Evolution:
    dρ_I/dt = -i[H₁,I(t), ρ_I(t)]
    where H₁,I(t) = U₀†(t) H₁(t) U₀(t)
"""

using DifferentialEquations
using Plots
using LinearAlgebra
using Random
using Distributions
using Statistics
using Interpolations

# ============================================================================
# Hamiltonian Decomposition
# ============================================================================

"""
For our system in the rotating frame:
    H_RF = (Δ + ζ(t))/2 * σz + f(t)/2 * σx
    
We split this as:
    H₀ = Δ/2 * σz              (free part - detuning)
    H₁(t) = ζ(t)/2 * σz + f(t)/2 * σx    (interaction part - noise + drive)
    
In the interaction picture, H₀ is transformed away and we only evolve
under the interaction-picture Hamiltonian H₁,I(t).
"""

# ============================================================================
# Interaction Picture Transformations
# ============================================================================

"""
    build_U0(t, Δ)

Build the free evolution operator U₀(t) = exp(-iH₀t).

For H₀ = (Δ/2)σz:
    U₀(t) = exp(-iΔt/2 * σz) = diag(exp(-iΔt/2), exp(iΔt/2))
"""
function build_U0(t, Δ)
    exp_factor_0 = exp(-im * Δ * t / 2)
    exp_factor_1 = exp(im * Δ * t / 2)
    return ComplexF64[exp_factor_0 0; 0 exp_factor_1]
end

"""
    transform_to_interaction_picture(H, t, Δ)

Transform Hamiltonian to interaction picture: H_I(t) = U₀†(t) H U₀(t)

For our system with H₁ = ζ(t)/2 * σz + f(t)/2 * σx:
    H₁,I(t) = U₀† H₁ U₀
"""
function transform_to_interaction_picture(H, t, Δ)
    U0 = build_U0(t, Δ)
    U0_dag = U0'
    return U0_dag * H * U0
end

"""
    build_interaction_hamiltonian(t, tg, θ0, Δ, ζ_t)

Build the interaction Hamiltonian H₁(t) in the Schrödinger picture.

H₁(t) = ζ(t)/2 * σz + f(t)/2 * σx
"""
function build_interaction_hamiltonian(t, tg, θ0, Δ, ζ_t)
    # Drive amplitude
    ft = (θ0/tg) * (1 - cos(2*π*t/tg))
    
    # Pauli matrices
    σz = ComplexF64[1 0; 0 -1]
    σx = ComplexF64[0 1; 1 0]
    
    # Interaction Hamiltonian (noise + drive)
    H1 = (ζ_t/2) * σz + (ft/2) * σx
    
    return H1
end

"""
    build_interaction_picture_hamiltonian(t, tg, θ0, Δ, ζ_t)

Build H₁,I(t) = U₀†(t) H₁(t) U₀(t) directly.

This is what actually drives the evolution in the interaction picture.
"""
function build_interaction_picture_hamiltonian(t, tg, θ0, Δ, ζ_t)
    # Build H₁ in Schrödinger picture
    H1 = build_interaction_hamiltonian(t, tg, θ0, Δ, ζ_t)
    
    # Transform to interaction picture
    H1_I = transform_to_interaction_picture(H1, t, Δ)
    
    return H1_I
end

# ============================================================================
# Vectorization for Superoperator Formalism
# ============================================================================

function vec_density_matrix(ρ::Matrix{<:Number})
    return [ρ[1,1], ρ[1,2], ρ[2,1], ρ[2,2]]
end

function unvec_density_matrix(ρ_vec::Vector{<:Number})
    return [ρ_vec[1] ρ_vec[2]; ρ_vec[3] ρ_vec[4]]
end

"""
    commutator_superoperator(H::Matrix)

Build superoperator for -i[H, ρ] acting on vectorized ρ.
"""
function commutator_superoperator(H::Matrix{<:Number})
    I2 = Matrix{ComplexF64}(I, 2, 2)
    L_comm = -im * (kron(I2, H) - kron(transpose(H), I2))
    return L_comm
end

# ============================================================================
# Noise Generation
# ============================================================================

function generate_ou_noise(τ_c, μ, σ, tspan, dt; u0=0.0, solver=LambaEulerHeun(), 
                          reltol=1e-6, abstol=1e-8)
    f(u, p, t) = (1/τ_c) * (μ - u)
    g(u, p, t) = sqrt(2*σ*σ/τ_c)
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, saveat=dt, abstol=abstol, reltol=reltol)
    
    return sol
end

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
# Averaged Liouvillian in Interaction Picture
# ============================================================================

"""
    compute_average_liouvillian_interaction_picture(times, tg, θ0, Δ, noise_ensemble)

Compute the time-dependent ensemble-averaged Liouvillian in the INTERACTION PICTURE.

This transforms out the detuning Δ, leaving only the averaged effect of the
drive + noise interaction.

# Arguments
- `times`: Time points
- `tg`: Gate time
- `θ0`: Rotation angle
- `Δ`: Detuning (qubit frequency - drive frequency)
- `noise_ensemble`: EnsembleSolution of noise trajectories

# Returns
- Vector of 4×4 averaged Liouvillian matrices ⟨L_I(t)⟩
"""
function compute_average_liouvillian_interaction_picture(times, tg, θ0, Δ, noise_ensemble)
    n_times = length(times)
    n_samples = length(noise_ensemble)
    
    # Store average Liouvillian at each time
    L_avg = [zeros(ComplexF64, 4, 4) for _ in 1:n_times]
    
    println("Computing averaged Liouvillian (Interaction Picture)...")
    println("  Detuning Δ = $Δ")
    println("  Averaging over $n_samples noise realizations...")
    
    for i in 1:n_samples
        # Get noise trajectory
        noise_sol = noise_ensemble[i]
        noise_interp = linear_interpolation(noise_sol.t, noise_sol[1,:], 
                                           extrapolation_bc=Flat())
        
        # Compute Liouvillian for this noise realization at each time
        for (j, t) in enumerate(times)
            ζ_t = noise_interp(t)
            
            # Build interaction-picture Hamiltonian
            H1_I = build_interaction_picture_hamiltonian(t, tg, θ0, Δ, ζ_t)
            
            # Convert to superoperator
            L_noisy = commutator_superoperator(H1_I)
            L_avg[j] += L_noisy
        end
        
        if i % 10 == 0
            println("    Processed $i/$n_samples realizations")
        end
    end
    
    # Average
    for j in 1:n_times
        L_avg[j] /= n_samples
    end
    
    println("  ✓ Averaged interaction-picture Liouvillian computed!")
    
    return L_avg
end

# ============================================================================
# Evolution with Interaction Picture Superoperator
# ============================================================================

"""
    evolve_interaction_picture(ρ0::Matrix, L_avg_I, times, Δ)

Evolve density matrix using averaged interaction-picture superoperator.

# Process:
1. Transform ρ₀ to interaction picture: ρ_I(0) = U₀†(0) ρ₀ U₀(0) = ρ₀
2. Evolve using ⟨L_I(t)⟩: dρ_I/dt = ⟨L_I(t)⟩ ρ_I
3. Transform back to Schrödinger picture: ρ_S(t) = U₀(t) ρ_I(t) U₀†(t)

# Arguments
- `ρ0`: Initial density matrix (Schrödinger picture)
- `L_avg_I`: Averaged interaction-picture Liouvillian from compute_average_liouvillian_interaction_picture
- `times`: Time points
- `Δ`: Detuning

# Returns
- Named tuple with ρ_evolved (Schrödinger picture), ρ_I_evolved (interaction picture), sol
"""
function evolve_interaction_picture(ρ0::Matrix, L_avg_I, times, Δ)
    n_times = length(times)
    
    # Initial state is same in both pictures at t=0
    ρ_I_0 = ComplexF64.(vec_density_matrix(ρ0))
    
    # ODE for interaction picture evolution
    function interaction_picture_ode!(dρ_I_vec, ρ_I_vec, p, t)
        times_ref, L_avg_I_ref = p
        # Find nearest time point
        t_idx = argmin(abs.(times_ref .- t))
        L_I_t = L_avg_I_ref[t_idx]
        
        # Interaction picture evolution
        dρ_I_vec .= L_I_t * ρ_I_vec
        return nothing
    end
    
    # Solve interaction picture evolution
    tspan = (times[1], times[end])
    prob = ODEProblem(interaction_picture_ode!, ρ_I_0, tspan, (times, L_avg_I))
    sol = solve(prob, Tsit5(), saveat=times)
    
    # Convert to density matrices (interaction picture)
    ρ_I_evolved = [unvec_density_matrix(sol.u[i]) for i in 1:length(sol.u)]
    
    # Transform back to Schrödinger picture: ρ_S(t) = U₀(t) ρ_I(t) U₀†(t)
    ρ_S_evolved = similar(ρ_I_evolved)
    for (i, t) in enumerate(times)
        U0 = build_U0(t, Δ)
        ρ_S_evolved[i] = U0 * ρ_I_evolved[i] * U0'
    end
    
    return (ρ_evolved=ρ_S_evolved, ρ_I_evolved=ρ_I_evolved, times=times, sol=sol, Δ=Δ)
end

# ============================================================================
# Main Workflow Functions
# ============================================================================

"""
    compute_averaged_superoperator_interaction_picture(tg, θ0, Δ, τ_c, σ, n_samples; n_times=1001)

Compute the averaged Liouvillian in the INTERACTION PICTURE.

This is the main function - compute it once, then apply to any initial state!

# Arguments
- `tg`: Gate time
- `θ0`: Rotation angle
- `Δ`: Detuning (ωq - ωd)
- `τ_c`: Noise correlation time
- `σ`: Noise amplitude
- `n_samples`: Number of noise realizations

# Returns
Named tuple with L_avg_I, times, Δ, noise_ensemble
"""
function compute_averaged_superoperator_interaction_picture(tg, θ0, Δ, τ_c, σ, n_samples; n_times=1001)
    println("\n" * "="^60)
    println("Interaction Picture: Computing Averaged Superoperator")
    println("="^60)
    println("Parameters:")
    println("  tg = $tg, θ0 = $θ0, Δ = $Δ")
    println("  τ_c = $τ_c, σ = $σ")
    println("  n_samples = $n_samples, n_times = $n_times")
    println()
    
    # Time grid
    times = collect(range(0, tg, length=n_times))
    dt = tg / 1000
    
    # Generate noise ensemble
    println("Step 1: Generating noise ensemble...")
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, (0.0, tg), dt, n_samples)
    println("  ✓ Noise ensemble generated!")
    println()
    
    # Compute averaged Liouvillian in interaction picture
    println("Step 2: Computing averaged interaction-picture Liouvillian...")
    L_avg_I = compute_average_liouvillian_interaction_picture(times, tg, θ0, Δ, noise_ensemble)
    println()
    
    println("✓ Interaction-picture superoperator ready!")
    println("  Use apply_superoperator_interaction_picture() to evolve states.")
    println()
    
    return (L_avg_I=L_avg_I, times=times, Δ=Δ, noise_ensemble=noise_ensemble,
            tg=tg, θ0=θ0, τ_c=τ_c, σ=σ)
end

"""
    apply_superoperator_interaction_picture(ρ0::Matrix, superop_I)

Apply the interaction-picture superoperator to evolve an initial density matrix.

# Arguments
- `ρ0`: Initial 2×2 density matrix
- `superop_I`: Result from compute_averaged_superoperator_interaction_picture()

# Returns
Named tuple with ρ_evolved, ρ_I_evolved, times, Δ
"""
function apply_superoperator_interaction_picture(ρ0::Matrix, superop_I)
    println("Evolving using interaction-picture superoperator...")
    result = evolve_interaction_picture(ρ0, superop_I.L_avg_I, superop_I.times, superop_I.Δ)
    println("  ✓ Evolution complete!")
    
    return result
end

# ============================================================================
# Helper Functions
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
    return ψ * ψ'
end

function extract_populations(ρ_evolved)
    P0 = [real(ρ[1,1]) for ρ in ρ_evolved]
    P1 = [real(ρ[2,2]) for ρ in ρ_evolved]
    return (P0=P0, P1=P1)
end

# ============================================================================
# Plotting
# ============================================================================

"""
    plot_interaction_picture_result(result; show_both=false)

Plot results from interaction-picture evolution.

If show_both=true, plots both Schrödinger and interaction picture populations.
"""
function plot_interaction_picture_result(result; show_both=false)
    pops_S = extract_populations(result.ρ_evolved)
    
    if show_both
        pops_I = extract_populations(result.ρ_I_evolved)
        
        p1 = plot(result.times, pops_S.P0, label="ρ₀₀ (Schrödinger)", 
                 xlabel="Time", ylabel="Population", linewidth=2,
                 title="Schrödinger Picture")
        plot!(p1, result.times, pops_S.P1, label="ρ₁₁", linewidth=2)
        
        p2 = plot(result.times, pops_I.P0, label="ρ₀₀ (Interaction)", 
                 xlabel="Time", ylabel="Population", linewidth=2,
                 title="Interaction Picture (Δ = $(result.Δ))")
        plot!(p2, result.times, pops_I.P1, label="ρ₁₁", linewidth=2)
        
        return plot(p1, p2, layout=(2,1), size=(800, 600))
    else
        p = plot(result.times, pops_S.P0, label="ρ₀₀ (Ground)", 
                xlabel="Time", ylabel="Population", linewidth=2,
                title="Interaction Picture Evolution (Δ = $(result.Δ))")
        plot!(p, result.times, pops_S.P1, label="ρ₁₁ (Excited)", linewidth=2)
        return p
    end
end

# ============================================================================
# Examples
# ============================================================================

"""
    example_interaction_picture()

Example: Evolution with interaction picture superoperator.
"""
function example_interaction_picture()
    println("\n" * "="^70)
    println("Example: Interaction Picture Formalism")
    println("="^70)
    
    # Parameters
    tg = 1.0
    θ0 = π/2
    Δ = 0.0  # On resonance (no detuning)
    τ_c = 0.1
    σ = 0.1
    n_samples = 50
    
    # STEP 1: Compute interaction-picture superoperator
    println("\n*** STEP 1: Computing superoperator (Interaction Picture) ***")
    sup_I = compute_averaged_superoperator_interaction_picture(tg, θ0, Δ, τ_c, σ, n_samples, n_times=501)
    
    # STEP 2: Apply to ground state
    println("\n*** STEP 2: Applying to ground state ***")
    result = apply_superoperator_interaction_picture(ground_state_matrix(), sup_I)
    
    # Final state
    ρ_final = result.ρ_evolved[end]
    println("\nFinal state (Schrödinger picture):")
    println("  ρ₀₀ = $(ρ_final[1,1])")
    println("  ρ₁₁ = $(ρ_final[2,2])")
    println()
    
    # Plot
    display(plot_interaction_picture_result(result, show_both=true))
    
    return sup_I
end

"""
    example_with_detuning()

Example: Show effect of detuning Δ ≠ 0.
"""
function example_with_detuning()
    println("\n" * "="^70)
    println("Example: Interaction Picture with Detuning")
    println("="^70)
    
    # Compute superoperators for different detunings
    Δ_values = [0.0, 0.5, 1.0]
    results = []
    
    for Δ in Δ_values
        println("\n--- Computing for Δ = $Δ ---")
        sup_I = compute_averaged_superoperator_interaction_picture(
            1.0, π/2, Δ, 0.1, 0.1, 30, n_times=501
        )
        result = apply_superoperator_interaction_picture(ground_state_matrix(), sup_I)
        push!(results, result)
    end
    
    # Plot comparison
    plots = []
    for (i, result) in enumerate(results)
        pops = extract_populations(result.ρ_evolved)
        p = plot(result.times, pops.P0, label="P₀", linewidth=2,
                title="Δ = $(Δ_values[i])", xlabel="Time", ylabel="Population")
        plot!(p, result.times, pops.P1, label="P₁", linewidth=2)
        push!(plots, p)
    end
    
    display(plot(plots..., layout=(3,1), size=(800, 900)))
    
    println("\n✓ Detuning changes the effective evolution in interaction picture!")
end

# ============================================================================
# Help Message
# ============================================================================

println("""
╔════════════════════════════════════════════════════════════╗
║  Interaction Picture (Dirac Picture) Module Loaded!       ║
╚════════════════════════════════════════════════════════════╝

KEY CONCEPT - Interaction Picture:
    Split H = H₀ + H₁ where:
    - H₀ = (Δ/2)σz  → Transforms away (free evolution)
    - H₁ = ζ(t)/2·σz + f(t)/2·σx  → Drives evolution

WORKFLOW:
    1. Compute superoperator (interaction picture):
       sup_I = compute_averaged_superoperator_interaction_picture(
           tg, θ0, Δ, τ_c, σ, n_samples
       )
    
    2. Apply to any initial state:
       result = apply_superoperator_interaction_picture(ρ0, sup_I)

ADVANTAGES:
    ✓ Transforms out detuning Δ
    ✓ Focuses on interaction dynamics
    ✓ Better for perturbative analysis
    ✓ Still compute once, apply to any ρ0!

Examples:
    example_interaction_picture()  - Basic usage
    example_with_detuning()        - Show Δ effect
""")