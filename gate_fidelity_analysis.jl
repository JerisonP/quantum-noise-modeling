"""
Gate Fidelity Analysis: Sweep Gate Time

Computes average gate fidelity as a function of gate time tg/τ_c
using the superoperator formalism with noise averaging.
"""

using DifferentialEquations
using Plots
using LinearAlgebra
using Random
using Distributions
using Statistics
using Interpolations

# Load the superoperator formalism
include("superoperator_formalism.jl")

# ============================================================================
# Pauli Eigenstate Definitions
# ============================================================================

"""
Six fiducial states (Pauli eigenstates):
- |z+⟩ = |0⟩ (eigenstate of σz with eigenvalue +1)
- |z-⟩ = |1⟩ (eigenstate of σz with eigenvalue -1)
- |x+⟩ = (|0⟩ + |1⟩)/√2 (eigenstate of σx with eigenvalue +1)
- |x-⟩ = (|0⟩ - |1⟩)/√2 (eigenstate of σx with eigenvalue -1)
- |y+⟩ = (|0⟩ + i|1⟩)/√2 (eigenstate of σy with eigenvalue +1)
- |y-⟩ = (|0⟩ - i|1⟩)/√2 (eigenstate of σy with eigenvalue -1)
"""

function pauli_eigenstates()
    # Define the six Pauli eigenstates as density matrices
    ρ_zp = ComplexF64[1 0; 0 0]  # |z+⟩⟨z+| = |0⟩⟨0|
    ρ_zm = ComplexF64[0 0; 0 1]  # |z-⟩⟨z-| = |1⟩⟨1|
    
    ρ_xp = ComplexF64[0.5 0.5; 0.5 0.5]  # |x+⟩⟨x+|
    ρ_xm = ComplexF64[0.5 -0.5; -0.5 0.5]  # |x-⟩⟨x-|
    
    ρ_yp = ComplexF64[0.5 -0.5im; 0.5im 0.5]  # |y+⟩⟨y+|
    ρ_ym = ComplexF64[0.5 0.5im; -0.5im 0.5]  # |y-⟩⟨y-|
    
    return [ρ_zp, ρ_zm, ρ_xp, ρ_xm, ρ_yp, ρ_ym]
end

# ============================================================================
# Ideal Gate (Unitary Evolution)
# ============================================================================

"""
    compute_ideal_gate(tg, θ0)

Compute the ideal unitary gate UI for a π/2 rotation without noise.

Returns the unitary matrix UI such that ρ_ideal = UI * ρ0 * UI†
"""
function compute_ideal_gate(tg, θ0)
    # For the shaped pulse f(t) = (θ0/tg)*(1-cos(2πt/tg))
    # The total rotation angle is: θ_total = ∫₀ᵗᵍ f(t)dt = θ0
    
    # For a rotation about x-axis by angle θ:
    # U = exp(-iθσx/2) = cos(θ/2)*I - i*sin(θ/2)*σx
    
    θ_total = θ0
    cos_half = cos(θ_total/2)
    sin_half = sin(θ_total/2)
    
    σx = ComplexF64[0 1; 1 0]
    UI = cos_half * I(2) - im * sin_half * σx
    
    return UI
end

# ============================================================================
# Average Fidelity Calculation
# ============================================================================

"""
    average_fidelity(ρ_ideal_list, ρ_noisy_list)

Compute average fidelity over 6 fiducial states.

F_avg = (1/6) Σₖ Tr(ρᵢᵈₑₐₗ₍ₖ₎ ρₙₒᵢₛy₍ₖ₎)

Returns: log₁₀(1 - F_avg) for plotting gate infidelity
"""
function average_fidelity(ρ_ideal_list, ρ_noisy_list)
    @assert length(ρ_ideal_list) == 6 "Must provide exactly 6 ideal states"
    @assert length(ρ_noisy_list) == 6 "Must provide exactly 6 noisy states"
    
    total_overlap = 0.0
    for k in 1:6
        # Compute Tr(ρᵢᵈₑₐₗ₍ₖ₎ * ρₙₒᵢₛy₍ₖ₎)
        overlap = tr(ρ_ideal_list[k] * ρ_noisy_list[k])
        total_overlap += real(overlap)
    end
    
    F_avg = total_overlap / 6
    
    # Return log₁₀(1 - F_avg) for infidelity
    infidelity = 1 - F_avg
    
    # Ensure infidelity is positive to avoid log of negative/zero
    if infidelity <= 0
        return -16.0  # Floor for numerical precision
    end
    
    return log10(infidelity)
end

# ============================================================================
# Gate Fidelity vs Gate Time
# ============================================================================

"""
    compute_gate_fidelity_vs_time(τ_c, θ0, σ, n_samples; tg_min_factor=0.1, tg_max_factor=10, n_points=20)

Sweep gate time tg and compute average gate fidelity.

# Arguments
- `τ_c`: Correlation time of noise
- `θ0`: Target rotation angle (typically π/2)
- `σ`: Noise amplitude
- `n_samples`: Number of noise realizations for averaging
- `tg_min_factor`: Minimum tg as fraction of τ_c (default: 0.1)
- `tg_max_factor`: Maximum tg as multiple of τ_c (default: 10)
- `n_points`: Number of points in the sweep

# Returns
Named tuple with (tg_values, tg_over_tau, log_infidelity)
"""
function compute_gate_fidelity_vs_time(τ_c, θ0, σ, n_samples; 
                                        tg_min_factor=0.1, 
                                        tg_max_factor=10, 
                                        n_points=20,
                                        n_times_per_gate=501,
                                        Δ=0.0)
    println("\n" * "="^70)
    println("Gate Fidelity vs Gate Time Analysis")
    println("="^70)
    println("Parameters:")
    println("  τ_c = $τ_c")
    println("  θ0 = $θ0 ($(θ0/π)π)")
    println("  σ = $σ")
    println("  Δ = $Δ")
    println("  n_samples = $n_samples noise realizations")
    println("  Sweeping tg from $(tg_min_factor)τ_c to $(tg_max_factor)τ_c")
    println("  Number of points: $n_points")
    println()
    
    # Generate tg values (log-spaced for better coverage)
    tg_values = 10 .^ range(log10(tg_min_factor * τ_c), 
                            log10(tg_max_factor * τ_c), 
                            length=n_points)
    tg_over_tau = tg_values ./ τ_c
    
    # Get Pauli eigenstates
    fiducial_states = pauli_eigenstates()
    
    # Storage for results
    log_infidelity = zeros(n_points)
    
    # Loop over gate times
    for (i, tg) in enumerate(tg_values)
        println("="^70)
        println("Point $i/$n_points: tg = $(round(tg, digits=4)) ($(round(tg/τ_c, digits=3))×τ_c)")
        println("="^70)
        
        # STEP 1: Compute ideal gate (no noise)
        UI = compute_ideal_gate(tg, θ0)
        ρ_ideal_list = [UI * ρ * UI' for ρ in fiducial_states]
        
        # STEP 2: Compute averaged superoperator with noise
        sup = compute_averaged_superoperator(tg, θ0, τ_c, σ, n_samples, 
                                             n_times=n_times_per_gate)
        
        # STEP 3: Apply to all 6 fiducial states
        println("\nApplying superoperator to 6 fiducial states...")
        ρ_noisy_list = []
        for (j, ρ0) in enumerate(fiducial_states)
            result = apply_superoperator(ρ0, sup)
            ρ_final = result.ρ_evolved[end]
            push!(ρ_noisy_list, ρ_final)
        end
        
        # STEP 4: Compute average fidelity
        log_infid = average_fidelity(ρ_ideal_list, ρ_noisy_list)
        log_infidelity[i] = log_infid
        
        println("\nResult: log₁₀(1 - F_avg) = $(round(log_infid, digits=4))")
        println()
    end
    
    println("="^70)
    println("Sweep Complete!")
    println("="^70)
    
    return (tg_values=tg_values, tg_over_tau=tg_over_tau, 
            log_infidelity=log_infidelity, τ_c=τ_c, σ=σ, θ0=θ0)
end

# ============================================================================
# Plotting
# ============================================================================

"""
    plot_fidelity_vs_time(result; kwargs...)

Plot log₁₀(1 - F_avg) vs tg/τ_c.
"""
function plot_fidelity_vs_time(result; 
                               title="Gate Fidelity vs Gate Time",
                               xlabel="Gate Time (tg/τc)",
                               ylabel="log₁₀(1 - F_avg)",
                               xscale=:log10,
                               linewidth=2,
                               markersize=4,
                               legend=:best)
    
    p = plot(result.tg_over_tau, result.log_infidelity,
            xlabel=xlabel, ylabel=ylabel, title=title,
            xscale=xscale, linewidth=linewidth,
            marker=:circle, markersize=markersize,
            label="σ = $(result.σ), θ0 = $(result.θ0/π)π",
            legend=legend, grid=true)
    
    return p
end

"""
    plot_fidelity_multiple_sigma(τ_c, θ0, σ_values, n_samples; kwargs...)

Compare fidelity for different noise amplitudes.
"""
function plot_fidelity_multiple_sigma(τ_c, θ0, σ_values, n_samples;
                                       tg_min_factor=0.1,
                                       tg_max_factor=10,
                                       n_points=20)
    println("\n" * "="^70)
    println("Comparing Multiple Noise Amplitudes")
    println("="^70)
    
    results = []
    
    for σ in σ_values
        println("\n\n*** Processing σ = $σ ***\n")
        result = compute_gate_fidelity_vs_time(τ_c, θ0, σ, n_samples,
                                               tg_min_factor=tg_min_factor,
                                               tg_max_factor=tg_max_factor,
                                               n_points=n_points)
        push!(results, result)
    end
    
    # Plot all on same axes
    p = plot(xlabel="Gate Time (tg/τc)", 
            ylabel="log₁₀(1 - F_avg)",
            title="Gate Fidelity vs Gate Time",
            xscale=:log10, linewidth=2,
            legend=:best, grid=true)
    
    for result in results
        plot!(p, result.tg_over_tau, result.log_infidelity,
             marker=:circle, markersize=3,
             label="σ = $(result.σ)")
    end
    
    display(p)
    
    return results
end

# ============================================================================
# Examples
# ============================================================================

"""
    example_fidelity_sweep()

Example: Sweep gate time and compute fidelity.
"""
function example_fidelity_sweep()
    # Parameters matching your code
    τ_c = 1.0
    θ0 = π/2
    σ = sqrt(1/(10000*τ_c))  # Your small noise case
    n_samples = 100  # Start with fewer for testing
    
    # Sweep from 0.1×τ_c to 10×τ_c
    result = compute_gate_fidelity_vs_time(τ_c, θ0, σ, n_samples,
                                           tg_min_factor=0.1,
                                           tg_max_factor=10,
                                           n_points=15)
    
    # Plot
    p = plot_fidelity_vs_time(result)
    display(p)
    
    return result
end

"""
    example_sigma_comparison()

Example: Compare different noise amplitudes.
"""
function example_sigma_comparison()
    τ_c = 1.0
    θ0 = π/2
    
    # Different noise amplitudes (matching your sweep)
    σ_values = [sqrt(1/(10000*τ_c)), sqrt(1/(1000*τ_c)), sqrt(1/(100*τ_c))]
    n_samples = 50
    
    results = plot_fidelity_multiple_sigma(τ_c, θ0, σ_values, n_samples,
                                          tg_min_factor=0.1,
                                          tg_max_factor=10,
                                          n_points=12)
    
    return results
end

"""
    example_full_sweep()

Full analysis: Match your original sweep σ from small to large.
"""
function example_full_sweep()
    τ_c = 1.0
    θ0 = π/2
    n_samples = 100
    
    # Generate 10 σ values from sqrt(1/(10000*τ_c)) to sqrt(1/τ_c)
    σ_min = sqrt(1/(10000*τ_c))
    σ_max = sqrt(1/τ_c)
    σ_values = 10 .^ range(log10(σ_min), log10(σ_max), length=10)
    
    # For each σ, compute fidelity at specific tg values
    tg_test = 0.7 * τ_c  # Your test case
    
    println("\n" * "="^70)
    println("Full Sweep: Multiple σ at fixed tg = $(tg_test)")
    println("="^70)
    
    log_infidelity_at_fixed_tg = []
    
    for σ in σ_values
        println("\n*** σ = $(round(σ, digits=6)) ***")
        
        # Compute for single tg value
        result = compute_gate_fidelity_vs_time(τ_c, θ0, σ, n_samples,
                                               tg_min_factor=tg_test/τ_c,
                                               tg_max_factor=tg_test/τ_c,
                                               n_points=1)
        
        push!(log_infidelity_at_fixed_tg, result.log_infidelity[1])
    end
    
    # Plot log_infidelity vs σ
    p = plot(σ_values, log_infidelity_at_fixed_tg,
            xlabel="Noise Amplitude σ", ylabel="log₁₀(1 - F_avg)",
            title="Gate Fidelity vs Noise (tg = $(tg_test))",
            xscale=:log10, linewidth=2, marker=:circle,
            markersize=4, legend=false, grid=true)
    
    display(p)
    
    return (σ_values=σ_values, log_infidelity=log_infidelity_at_fixed_tg)
end

# ============================================================================
# Help Message
# ============================================================================

println("""
╔════════════════════════════════════════════════════════════╗
║  Gate Fidelity Analysis Module Loaded!                    ║
╚════════════════════════════════════════════════════════════╝

MAIN WORKFLOW:

1. Sweep gate time tg from 0.1×τ_c to 10×τ_c:
   result = compute_gate_fidelity_vs_time(τ_c, θ0, σ, n_samples)
   plot_fidelity_vs_time(result)

2. Compare multiple noise amplitudes:
   results = plot_fidelity_multiple_sigma(τ_c, θ0, [σ1, σ2, σ3], n_samples)

3. Full sweep over σ at fixed tg:
   result = example_full_sweep()

EXAMPLES:
   example_fidelity_sweep()    - Single σ, sweep tg
   example_sigma_comparison()  - Multiple σ, sweep tg
   example_full_sweep()        - Sweep σ at fixed tg

KEY FUNCTIONS:
   compute_gate_fidelity_vs_time()  - Main analysis
   average_fidelity()               - Fidelity calculation
   pauli_eigenstates()              - Get 6 fiducial states
   compute_ideal_gate()             - Ideal unitary
""")