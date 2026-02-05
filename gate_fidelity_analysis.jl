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

F_avg = (1/d²) Σₖ Tr(ρᵢᵈₑₐₗ₍ₖ₎ · ρₙₒᵢₛy₍ₖ₎ · ρᵢᵈₑₐₗ₍ₖ₎†)

For qubits, d=2, so we use d²=4 but average over 6 Pauli eigenstates:
F_avg = (1/6) Σₖ Tr(ρᵢᵈₑₐₗ₍ₖ₎ · ρₙₒᵢₛy₍ₖ₎ · ρᵢᵈₑₐₗ₍ₖ₎†)

Returns: log₁₀(1 - F_avg) for plotting gate infidelity
"""
function average_fidelity(ρ_ideal_list, ρ_noisy_list; debug=false)
    @assert length(ρ_ideal_list) == 6 "Must provide exactly 6 ideal states"
    @assert length(ρ_noisy_list) == 6 "Must provide exactly 6 noisy states"
    
    total_overlap = 0.0
    
    if debug
        println("\nDEBUG: Fidelity calculation details:")
        state_names = ["z+", "z-", "x+", "x-", "y+", "y-"]
    end
    
    for k in 1:6
        # Compute Tr(ρᵢᵈₑₐₗ₍ₖ₎ · ρₙₒᵢₛy₍ₖ₎ · ρᵢᵈₑₐₗ₍ₖ₎†)
        ρ_ideal_dag = ρ_ideal_list[k]'
        overlap = tr(ρ_ideal_list[k] * ρ_noisy_list[k] * ρ_ideal_dag)
        overlap_real = real(overlap)
        total_overlap += overlap_real
        
        if debug
            println("  State $(state_names[k]):")
            println("    Tr(ρ_ideal · ρ_noisy · ρ_ideal†) = $(round(overlap, digits=6))")
            println("    Real part = $(round(overlap_real, digits=6))")
            println("    ρ_ideal trace = $(round(tr(ρ_ideal_list[k]), digits=6))")
            println("    ρ_noisy trace = $(round(tr(ρ_noisy_list[k]), digits=6))")
            println("    ρ_ideal is Hermitian? $(ishermitian(ρ_ideal_list[k]))")
            println("    ρ_noisy is Hermitian? $(ishermitian(ρ_noisy_list[k]))")
        end
    end
    
    F_avg = total_overlap / 6
    
    if debug
        println("\n  Total overlap = $(round(total_overlap, digits=6))")
        println("  F_avg = $(round(F_avg, digits=6))")
        println("  1 - F_avg = $(round(1 - F_avg, digits=6))")
    end
    
    # Return log₁₀(1 - F_avg) for infidelity
    infidelity = 1 - F_avg
    
    # Ensure infidelity is positive to avoid log of negative/zero
    if infidelity <= 0
        if debug
            println("  WARNING: infidelity <= 0, using floor value -16")
        end
        return -16.0  # Floor for numerical precision
    end
    
    log_infid = log10(infidelity)
    
    if debug
        println("  log₁₀(1 - F_avg) = $(round(log_infid, digits=6))")
    end
    
    return log_infid
end

# ============================================================================
# Gate Fidelity vs Noise Amplitude (Sweeping σ)
# ============================================================================

"""
    compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples; σ_min_factor=1/10000, σ_max_factor=1, n_points=20)

Sweep noise amplitude σ and compute average gate fidelity at FIXED gate time tg.

# Arguments
- `τ_c`: Correlation time of noise
- `tg`: Gate time (FIXED - typically 0.7×τ_c)
- `θ0`: Target rotation angle (typically π/2)
- `n_samples`: Number of noise realizations for averaging
- `σ_min_factor`: Minimum σ² as 1/(factor×τ_c), default: 1/10000
- `σ_max_factor`: Maximum σ² as 1/(factor×τ_c), default: 1
- `n_points`: Number of points in the sweep

# Returns
Named tuple with (σ_values, log_infidelity, τ_c, tg, θ0)
"""
function compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples; 
                                        σ_min_factor=10000,  # σ² = 1/(10000×τ_c)
                                        σ_max_factor=1,      # σ² = 1/τ_c
                                        n_points=20,
                                        n_times_per_gate=501,
                                        Δ=0.0)
    println("\n" * "="^70)
    println("Gate Fidelity vs Noise Amplitude Analysis")
    println("="^70)
    println("Parameters:")
    println("  τ_c = $τ_c")
    println("  tg = $tg ($(round(tg/τ_c, digits=3))×τ_c) [FIXED]")
    println("  θ0 = $θ0 ($(θ0/π)π)")
    println("  Δ = $Δ")
    println("  n_samples = $n_samples noise realizations")
    println("  Sweeping σ from √(1/$(σ_min_factor)τ_c) to √(1/$(σ_max_factor)τ_c)")
    println("  Number of points: $n_points")
    println()
    
    # Generate σ values (log-spaced)
    # σ² ranges from 1/(σ_min_factor×τ_c) to 1/(σ_max_factor×τ_c)
    σ_min = sqrt(1/(σ_min_factor * τ_c))
    σ_max = sqrt(1/(σ_max_factor * τ_c))
    
    σ_values = 10 .^ range(log10(σ_min), log10(σ_max), length=n_points)
    
    # Get Pauli eigenstates
    fiducial_states = pauli_eigenstates()
    
    # STEP 1: Compute ideal gate ONCE (no noise, same for all σ)
    println("Computing ideal gate (no noise)...")
    UI = compute_ideal_gate(tg, θ0)
    ρ_ideal_list = [UI * ρ * UI' for ρ in fiducial_states]
    println("✓ Ideal gate computed")
    println()
    
    # Storage for results
    log_infidelity = zeros(n_points)
    
    # Loop over noise amplitudes
    for (i, σ) in enumerate(σ_values)
        println("="^70)
        println("Point $i/$n_points: σ = $(round(σ, sigdigits=4))")
        println("  (σ² = $(round(σ^2, sigdigits=4)), σ²×τ_c = $(round(σ^2*τ_c, sigdigits=4)))")
        println("="^70)
        
        # STEP 2: Compute averaged superoperator with this noise amplitude
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
        debug_mode = (i == 1)  # Debug first point
        log_infid = average_fidelity(ρ_ideal_list, ρ_noisy_list, debug=debug_mode)
        log_infidelity[i] = log_infid
        
        # Also compute the actual fidelity for diagnostics
        F_avg = 10^(-log_infid) > 1.0 ? 0.0 : 1.0 - 10^(log_infid)
        
        if !debug_mode
            println("\nResult:")
            println("  F_avg = $(round(F_avg, digits=6))")
            println("  1 - F_avg = $(round(1-F_avg, digits=6))")
            println("  log₁₀(1 - F_avg) = $(round(log_infid, digits=4))")
            println("  σ²×τ_c = $(round(σ^2 * τ_c, sigdigits=4))")
        end
        println()
    end
    
    println("="^70)
    println("Sweep Complete!")
    println("="^70)
    
    return (σ_values=σ_values, log_infidelity=log_infidelity, 
            τ_c=τ_c, tg=tg, θ0=θ0)
end

# ============================================================================
# Plotting
# ============================================================================

"""
    plot_fidelity_vs_sigma(result; kwargs...)

Plot log₁₀(1 - F_avg) vs dimensionless noise parameter σ²×τ_c.
"""
function plot_fidelity_vs_sigma(result; 
                                title="Gate Fidelity vs Noise Amplitude",
                                xlabel="Noise Parameter (σ²×τc)",
                                ylabel="log₁₀(1 - F_avg)",
                                xscale=:log10,
                                linewidth=2,
                                markersize=4,
                                legend=:best)
    
    # Compute dimensionless noise parameter
    noise_param = (result.σ_values .^ 2) .* result.τ_c
    
    p = plot(noise_param, result.log_infidelity,
            xlabel=xlabel, ylabel=ylabel, title=title,
            xscale=xscale, linewidth=linewidth,
            marker=:circle, markersize=markersize,
            label="tg = $(result.tg) ($(round(result.tg/result.τ_c, digits=2))×τc)",
            legend=legend, grid=true)
    
    return p
end

"""
    plot_fidelity_multiple_tg(τ_c, θ0, tg_values, n_samples; kwargs...)

Compare fidelity for different gate times (sweeping σ for each).
"""
function plot_fidelity_multiple_tg(τ_c, θ0, tg_values, n_samples;
                                   σ_min_factor=10000,
                                   σ_max_factor=1,
                                   n_points=20)
    println("\n" * "="^70)
    println("Comparing Multiple Gate Times")
    println("="^70)
    
    results = []
    
    for tg in tg_values
        println("\n\n*** Processing tg = $tg ($(round(tg/τ_c, digits=2))×τ_c) ***\n")
        result = compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples,
                                                σ_min_factor=σ_min_factor,
                                                σ_max_factor=σ_max_factor,
                                                n_points=n_points)
        push!(results, result)
    end
    
    # Plot all on same axes
    p = plot(xlabel="Noise Parameter (σ²×τc)", 
            ylabel="log₁₀(1 - F_avg)",
            title="Gate Fidelity vs Noise Amplitude",
            xscale=:log10, linewidth=2,
            legend=:best, grid=true)
    
    for result in results
        noise_param = (result.σ_values .^ 2) .* result.τ_c
        plot!(p, noise_param, result.log_infidelity,
             marker=:circle, markersize=3,
             label="tg = $(round(result.tg/result.τ_c, digits=2))×τc")
    end
    
    display(p)
    
    return results
end

# ============================================================================
# Examples
# ============================================================================

"""
    example_fidelity_sweep()

Example: Sweep noise amplitude σ at fixed gate time (matching your original code).
"""
function example_fidelity_sweep()
    # Parameters matching your code
    τ_c = 1.0
    tg = 0.7 * τ_c  # Fixed gate time
    θ0 = π/2
    n_samples = 100  # Start with fewer for testing
    
    # Sweep σ from sqrt(1/(10000×τ_c)) to sqrt(1/τ_c)
    result = compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples,
                                            σ_min_factor=10000,
                                            σ_max_factor=1,
                                            n_points=15)
    
    # Plot
    p = plot_fidelity_vs_sigma(result)
    display(p)
    
    return result
end

"""
    example_tg_comparison()

Example: Compare different gate times (sweeping σ for each).
"""
function example_tg_comparison()
    τ_c = 1.0
    θ0 = π/2
    
    # Different gate times
    tg_values = [0.1*τ_c, 0.7*τ_c, 2.0*τ_c]
    n_samples = 50
    
    results = plot_fidelity_multiple_tg(τ_c, θ0, tg_values, n_samples,
                                       σ_min_factor=10000,
                                       σ_max_factor=1,
                                       n_points=12)
    
    return results
end

"""
    example_full_sweep()

Full analysis matching your original code exactly.
Sweeps σ from sqrt(1/(10000×τ_c)) to sqrt(1/τ_c) at tg = 0.7×τ_c
"""
function example_full_sweep()
    τ_c = 1.0
    tg = 0.7 * τ_c  # Your test case
    θ0 = π/2
    n_samples = 1000  # Your original: 1000 samples
    
    println("\n" * "="^70)
    println("Full Sweep: Matching Original Code")
    println("  τ_c = $τ_c")
    println("  tg = $tg ($(tg/τ_c)×τ_c)")
    println("  θ0 = $θ0")
    println("  Sweeping σ from √(1/10000τ_c) to √(1/τ_c)")
    println("  Using $n_samples noise realizations per point")
    println("="^70)
    
    # Sweep with 100 points (your original had 100)
    result = compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples,
                                            σ_min_factor=10000,
                                            σ_max_factor=1,
                                            n_points=100)
    
    # Plot
    p = plot_fidelity_vs_sigma(result,
                               title="Gate Fidelity vs Noise (Full Sweep)")
    display(p)
    
    return result
end

"""
    test_ideal_gate()

Test that the ideal gate gives perfect fidelity (F_avg = 1.0) with itself.
"""
function test_ideal_gate()
    println("\n" * "="^70)
    println("Testing Ideal Gate (No Noise)")
    println("="^70)
    
    τ_c = 1.0
    tg = 0.1 * τ_c
    θ0 = π/2
    
    # Get Pauli eigenstates
    fiducial_states = pauli_eigenstates()
    
    # Compute ideal gate
    UI = compute_ideal_gate(tg, θ0)
    println("\nIdeal gate UI:")
    display(UI)
    println()
    
    # Apply to fiducial states
    ρ_ideal_list = [UI * ρ * UI' for ρ in fiducial_states]
    
    # Use same states for "noisy" (should give F_avg = 1.0)
    ρ_noisy_list = ρ_ideal_list
    
    # Compute fidelity
    println("\nComputing fidelity with itself (should be 1.0):")
    log_infid = average_fidelity(ρ_ideal_list, ρ_noisy_list, debug=true)
    
    F_avg = 1.0 - 10^(log_infid)
    
    println("\n" * "="^70)
    println("Result:")
    println("  F_avg = $(round(F_avg, digits=10))")
    println("  Expected: 1.0")
    println("  Match? $(isapprox(F_avg, 1.0, atol=1e-10))")
    println("="^70)
    
    return F_avg
end

"""
    run_specific_tg(tg_factor; n_samples=100, n_points=20)

Quick function to run analysis at a specific tg = tg_factor × τ_c.

# Examples
```julia
# tg = 0.1×τ_c (fast gate)
result = run_specific_tg(0.1)

# tg = 0.7×τ_c (your original)
result = run_specific_tg(0.7, n_samples=1000, n_points=100)

# tg = 10×τ_c (slow gate)
result = run_specific_tg(10.0)
```
"""
function run_specific_tg(tg_factor; τ_c=1.0, θ0=π/2, n_samples=100, n_points=20)
    tg = tg_factor * τ_c
    
    println("\n" * "="^70)
    println("Gate Fidelity Analysis: tg = $(tg_factor)×τ_c")
    println("="^70)
    
    result = compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples,
                                            σ_min_factor=10000,
                                            σ_max_factor=1,
                                            n_points=n_points)
    
    # Plot
    p = plot_fidelity_vs_sigma(result,
                               title="Gate Fidelity: tg = $(tg_factor)×τc")
    display(p)
    
    return result
end

# ============================================================================
# Help Message
# ============================================================================

println("""
╔════════════════════════════════════════════════════════════╗
║  Gate Fidelity Analysis Module Loaded!                    ║
╚════════════════════════════════════════════════════════════╝

MAIN WORKFLOW (Matching Your Original Code):

1. Sweep σ from √(1/10000τ_c) to √(1/τ_c) at FIXED tg:
   result = compute_gate_fidelity_vs_sigma(τ_c, tg, θ0, n_samples)
   plot_fidelity_vs_sigma(result)

2. Quick analysis at specific tg:
   result = run_specific_tg(0.1)    # tg = 0.1×τ_c
   result = run_specific_tg(0.7)    # tg = 0.7×τ_c
   result = run_specific_tg(10.0)   # tg = 10×τ_c

3. Compare multiple gate times:
   results = plot_fidelity_multiple_tg(τ_c, θ0, [0.1, 0.7, 2.0], n_samples)

EXAMPLES:
   run_specific_tg(0.1)        - Fast gate (tg = 1/10 τ_c)
   run_specific_tg(0.7)        - Medium gate (tg = 0.7 τ_c)
   example_fidelity_sweep()    - Default sweep
   example_tg_comparison()     - Compare tg values
   example_full_sweep()        - Full 100-point sweep

TYPICAL USAGE:
   # Quick test with tg = 1/10 τ_c
   result = run_specific_tg(0.1, n_samples=100, n_points=15)
   
   # Production run with tg = 1/10 τ_c
   result = run_specific_tg(0.1, n_samples=1000, n_points=100)

KEY FUNCTIONS:
   compute_gate_fidelity_vs_sigma()  - Main analysis (sweeps σ)
   run_specific_tg()                 - Convenience wrapper
   average_fidelity()                - Fidelity calculation
""")