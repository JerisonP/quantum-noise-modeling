"""
Density Matrix Evolution with Ornstein-Uhlenbeck Noise

Combines the working density matrix simulation with OU noise modeling.
"""

using DifferentialEquations
using Plots
using LinearAlgebra
using Random
using Distributions
using Statistics
using Interpolations

# ============================================================================
# Core ODE Function (Without Noise)
# ============================================================================

"""
    density_matrix_ode!(du, u, p, t)

Time evolution of density matrix with time-dependent drive (no noise).
"""
function density_matrix_ode!(du, u, p, t)
    # Unpack the components of the density matrix
    ρ00Re, ρ00Im, ρ01Re, ρ01Im, ρ10Re, ρ10Im, ρ11Re, ρ11Im = u
    
    # Unpack parameters
    tg, θ0 = p
    
    # Calculate the time-dependent drive amplitude
    ft = (θ0/tg) * (1 - cos((2*π*t)/tg))
    
    # Define the real matrix with previously provided simplifications
    matRe = [
        [0.5*(-ρ10Im + ρ01Im)*ft, 0.5*(ρ00Im - ρ11Im)*ft],
        [0.5*(ρ11Im - ρ00Im)*ft, 0.5*(ρ10Im - ρ01Im)*ft]
    ]
    
    # Define the imaginary matrix based on the derived terms
    matIm = [
        [0.5*(ρ10Re - ρ01Re)*ft, 0.5*(ρ11Re - ρ00Re)*ft],
        [0.5*(ρ00Re - ρ11Re)*ft, 0.5*(ρ01Re - ρ10Re)*ft]
    ]
    
    # Assign derivatives based on matRe and matIm entries
    du[1] = matRe[1][1]  # ρ̇00Re
    du[2] = matIm[1][1]  # ρ̇00Im
    du[3] = matRe[1][2]  # ρ̇01Re
    du[4] = matIm[1][2]  # ρ̇01Im
    du[5] = matRe[2][1]  # ρ̇10Re
    du[6] = matIm[2][1]  # ρ̇10Im
    du[7] = matRe[2][2]  # ρ̇11Re
    du[8] = matIm[2][2]  # ρ̇11Im
    
    return nothing
end

# ============================================================================
# Core ODE Function (With Noise)
# ============================================================================

"""
    density_matrix_ode_with_noise!(du, u, p, t)

Time evolution of density matrix with time-dependent drive AND noise trajectory.

# Parameters p:
- p[1]: tg (gate time)
- p[2]: θ0 (rotation angle)
- p[3]: noise_trajectory (interpolated noise function)
"""
function density_matrix_ode_with_noise!(du, u, p, t)
    # Unpack the components of the density matrix
    ρ00Re, ρ00Im, ρ01Re, ρ01Im, ρ10Re, ρ10Im, ρ11Re, ρ11Im = u
    
    # Unpack parameters
    tg, θ0, noise_interp = p
    
    # Calculate the time-dependent drive amplitude
    ft = (θ0/tg) * (1 - cos((2*π*t)/tg))
    
    # Get noise value at current time
    ζ_t = noise_interp(t)
    
    # Modified drive with noise
    ft_noisy = ft * (1 + ζ_t)
    
    # Define the real matrix with noise-modified drive
    matRe = [
        [0.5*(-ρ10Im + ρ01Im)*ft_noisy, 0.5*(ρ00Im - ρ11Im)*ft_noisy],
        [0.5*(ρ11Im - ρ00Im)*ft_noisy, 0.5*(ρ10Im - ρ01Im)*ft_noisy]
    ]
    
    # Define the imaginary matrix with noise-modified drive
    matIm = [
        [0.5*(ρ10Re - ρ01Re)*ft_noisy, 0.5*(ρ11Re - ρ00Re)*ft_noisy],
        [0.5*(ρ00Re - ρ11Re)*ft_noisy, 0.5*(ρ01Re - ρ10Re)*ft_noisy]
    ]
    
    # Assign derivatives
    du[1] = matRe[1][1]
    du[2] = matIm[1][1]
    du[3] = matRe[1][2]
    du[4] = matIm[1][2]
    du[5] = matRe[2][1]
    du[6] = matIm[2][1]
    du[7] = matRe[2][2]
    du[8] = matIm[2][2]
    
    return nothing
end

# ============================================================================
# Ornstein-Uhlenbeck Noise Generation
# ============================================================================

"""
    generate_ou_noise(τ_c, μ, σ, tspan, dt; u0=0.0, solver=LambaEulerHeun())

Generate a single Ornstein-Uhlenbeck noise trajectory.

# Arguments
- `τ_c`: Correlation time
- `μ`: Mean value
- `σ`: Standard deviation
- `tspan`: Time span (start, end)
- `dt`: Time step for saving
- `u0`: Initial value (default: 0.0)
- `solver`: SDE solver (default: LambaEulerHeun())

# Returns
- Solution object with noise trajectory
"""
function generate_ou_noise(τ_c, μ, σ, tspan, dt; u0=0.0, solver=LambaEulerHeun(), 
                          reltol=1e-6, abstol=1e-8)
    # Drift term: mean reversion
    f(u, p, t) = (1/τ_c) * (μ - u)
    
    # Diffusion term: stochastic component
    g(u, p, t) = sqrt(2*σ*σ/τ_c)
    
    # Create and solve SDE
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, saveat=dt, abstol=abstol, reltol=reltol)
    
    return sol
end

"""
    generate_ou_ensemble(τ_c, μ, σ, tspan, dt, n_samples; kwargs...)

Generate an ensemble of Ornstein-Uhlenbeck noise trajectories.

# Arguments
- `τ_c`: Correlation time
- `μ`: Mean value  
- `σ`: Standard deviation
- `tspan`: Time span
- `dt`: Time step
- `n_samples`: Number of trajectories to generate

# Returns
- EnsembleSolution with all trajectories
"""
function generate_ou_ensemble(τ_c, μ, σ, tspan, dt, n_samples; 
                              solver=LambaEulerHeun(), reltol=1e-6, abstol=1e-8)
    # Drift and diffusion functions
    f(u, p, t) = (1/τ_c) * (μ - u)
    g(u, p, t) = sqrt(2*σ*σ/τ_c)
    
    # Problem function: random initial conditions
    function prob_func(prob, i, repeat)
        remake(prob, u0=rand(Normal(0, σ)))
    end
    
    # Create ensemble problem
    prob = SDEProblem(f, g, 0.0, tspan)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    
    # Solve ensemble
    sim = solve(ensemble_prob, solver, EnsembleThreads(), 
                trajectories=n_samples, saveat=dt, abstol=abstol, reltol=reltol)
    
    return sim
end

"""
    autocorrelation_function(trajectories_matrix, mean_trajectory, max_lag)

Compute normalized autocorrelation function of noise ensemble.
"""
function autocorrelation_function(trajectories_matrix, mean_trajectory, max_lag)
    n = size(trajectories_matrix, 1)  # number of time points
    covariances = Vector{Float64}(undef, max_lag+1)
    diff_matrix = trajectories_matrix .- mean_trajectory
    
    for τ in 0:max_lag
        cov_τ = 0.0
        for t in 1:(n - τ)
            cov_τ += mean(diff_matrix[t, :] .* diff_matrix[t+τ, :])
        end
        covariances[τ + 1] = cov_τ / (n - τ)
    end
    
    return covariances / covariances[1]
end

# ============================================================================
# Initial State Functions
# ============================================================================

function ground_state()
    return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end

function excited_state()
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
end

function superposition_state(α::Number, β::Number)
    norm = sqrt(abs2(α) + abs2(β))
    α, β = α/norm, β/norm
    
    ρ00 = abs2(α)
    ρ01 = α * conj(β)
    ρ10 = conj(α) * β
    ρ11 = abs2(β)
    
    return [real(ρ00), imag(ρ00), real(ρ01), imag(ρ01),
            real(ρ10), imag(ρ10), real(ρ11), imag(ρ11)]
end

# ============================================================================
# Simulation Functions (No Noise)
# ============================================================================

"""
    simulate(u0, tg, θ0; solver=Tsit5(), reltol=1e-10, abstol=1e-12)

Simulate density matrix evolution WITHOUT noise.
"""
function simulate(u0, tg, θ0; solver=Tsit5(), reltol=1e-10, abstol=1e-12)
    p = [tg, θ0]
    tspan = (0.0, tg)
    prob = ODEProblem(density_matrix_ode!, u0, tspan, p)
    sol = solve(prob, solver, reltol=reltol, abstol=abstol)
    return sol
end

# ============================================================================
# Simulation Functions (With Noise)
# ============================================================================

"""
    simulate_with_noise(u0, tg, θ0, noise_sol; solver=Tsit5(), reltol=1e-10, abstol=1e-12)

Simulate density matrix evolution WITH a given noise trajectory.

# Arguments
- `u0`: Initial state
- `tg`: Gate time
- `θ0`: Rotation angle
- `noise_sol`: Noise trajectory solution (from generate_ou_noise)
- `solver`: ODE solver
- `reltol`, `abstol`: Tolerances

# Returns
- Solution object
"""
function simulate_with_noise(u0, tg, θ0, noise_sol; 
                            solver=Tsit5(), reltol=1e-10, abstol=1e-12)
    # Create interpolation function for noise
    noise_interp = linear_interpolation(noise_sol.t, noise_sol[1,:], extrapolation_bc=Flat())
    
    # Pack parameters with noise interpolation (use tuple for better performance)
    p = (tg, θ0, noise_interp)
    tspan = (0.0, tg)
    
    # Solve with noise
    prob = ODEProblem(density_matrix_ode_with_noise!, u0, tspan, p)
    sol = solve(prob, solver, reltol=reltol, abstol=abstol)
    
    return sol
end

"""
    simulate_noisy_ensemble(u0, tg, θ0, τ_c, σ, n_samples; dt=nothing, kwargs...)

Simulate an ensemble of noisy density matrix evolutions.

# Arguments
- `u0`: Initial state
- `tg`: Gate time
- `θ0`: Rotation angle
- `τ_c`: Noise correlation time
- `σ`: Noise standard deviation
- `n_samples`: Number of ensemble members
- `dt`: Time step (default: tg/1000)

# Returns
- Array of solution objects, one for each noise realization
"""
function simulate_noisy_ensemble(u0, tg, θ0, τ_c, σ, n_samples; 
                                dt=nothing, solver=Tsit5(), 
                                reltol=1e-10, abstol=1e-12)
    if dt === nothing
        dt = tg / 1000
    end
    
    tspan = (0.0, tg)
    
    # Generate noise ensemble
    noise_ensemble = generate_ou_ensemble(τ_c, 0.0, σ, tspan, dt, n_samples)
    
    # Simulate density matrix for each noise trajectory
    solutions = []
    
    for i in 1:n_samples
        noise_sol = noise_ensemble[i]
        sol = simulate_with_noise(u0, tg, θ0, noise_sol, 
                                 solver=solver, reltol=reltol, abstol=abstol)
        push!(solutions, sol)
    end
    
    return solutions
end

# ============================================================================
# Analysis Functions
# ============================================================================

function extract_density_matrix(sol)
    ρ00 = sol[1,:] .+ im .* sol[2,:]
    ρ01 = sol[3,:] .+ im .* sol[4,:]
    ρ10 = sol[5,:] .+ im .* sol[6,:]
    ρ11 = sol[7,:] .+ im .* sol[8,:]
    return (ρ00=ρ00, ρ01=ρ01, ρ10=ρ10, ρ11=ρ11)
end

function get_populations(sol)
    P0 = sol[1,:]
    P1 = sol[7,:]
    return (P0=P0, P1=P1)
end

function check_trace(sol)
    return sol[1,:] .+ sol[7,:]
end

function final_state(sol)
    return [sol[1,end]+im*sol[2,end]  sol[3,end]+im*sol[4,end];
            sol[5,end]+im*sol[6,end]  sol[7,end]+im*sol[8,end]]
end

"""
    ensemble_statistics(solutions)

Compute mean and std of populations across ensemble.

# Returns
Named tuple with (P0_mean, P0_std, P1_mean, P1_std, times)
"""
function ensemble_statistics(solutions)
    n_samples = length(solutions)
    
    # Find common time grid (use the first solution's times as reference)
    times = solutions[1].t
    n_times = length(times)
    
    # Collect all populations (interpolate to common time grid)
    P0_all = zeros(n_times, n_samples)
    P1_all = zeros(n_times, n_samples)
    
    for (i, sol) in enumerate(solutions)
        # If solution has different time points, interpolate
        if length(sol.t) != n_times || sol.t != times
            # Create interpolations
            P0_interp = linear_interpolation(sol.t, sol[1,:], extrapolation_bc=Flat())
            P1_interp = linear_interpolation(sol.t, sol[7,:], extrapolation_bc=Flat())
            
            # Evaluate at common times
            P0_all[:, i] = [P0_interp(t) for t in times]
            P1_all[:, i] = [P1_interp(t) for t in times]
        else
            # Direct assignment if times match
            P0_all[:, i] = sol[1,:]
            P1_all[:, i] = sol[7,:]
        end
    end
    
    # Compute statistics
    P0_mean = mean(P0_all, dims=2)[:]
    P0_std = std(P0_all, dims=2)[:]
    P1_mean = mean(P1_all, dims=2)[:]
    P1_std = std(P1_all, dims=2)[:]
    
    return (P0_mean=P0_mean, P0_std=P0_std, 
            P1_mean=P1_mean, P1_std=P1_std, times=times)
end

# ============================================================================
# Plotting Functions
# ============================================================================

function plot_populations(sol; title="", legend=:best)
    p = plot(sol.t, sol[1,:], label="ρ₀₀ (Ground State)",
             xlabel="Time", ylabel="Population",
             linewidth=2, title=title, legend=legend)
    plot!(p, sol.t, sol[7,:], label="ρ₁₁ (Excited State)", linewidth=2)
    return p
end

function plot_coherences(sol; title="", legend=:best)
    p = plot(sol.t, sol[3,:], label="Re(ρ₀₁)",
             xlabel="Time", ylabel="Coherence",
             linewidth=2, title=title, legend=legend)
    plot!(p, sol.t, sol[4,:], label="Im(ρ₀₁)", linewidth=2)
    return p
end

"""
    plot_ensemble(solutions; n_plot=10, title="Noisy Ensemble")

Plot population dynamics for multiple ensemble members.

# Arguments
- `solutions`: Array of solution objects
- `n_plot`: Number of trajectories to plot (default: 10)
- `title`: Plot title
"""
function plot_ensemble(solutions; n_plot=10, title="Noisy Ensemble")
    n_plot = min(n_plot, length(solutions))
    
    p = plot(xlabel="Time", ylabel="Population", 
             title=title, legend=:outertopright)
    
    # Plot subset of trajectories
    for i in 1:n_plot
        plot!(p, solutions[i].t, solutions[i][1,:], 
              label=(i==1 ? "ρ₀₀" : ""), alpha=0.3, color=:blue, linewidth=1)
        plot!(p, solutions[i].t, solutions[i][7,:], 
              label=(i==1 ? "ρ₁₁" : ""), alpha=0.3, color=:red, linewidth=1)
    end
    
    return p
end

"""
    plot_ensemble_with_stats(solutions; title="Ensemble with Statistics")

Plot ensemble with mean and standard deviation bands.
"""
function plot_ensemble_with_stats(solutions; title="Ensemble with Statistics")
    stats = ensemble_statistics(solutions)
    
    p = plot(xlabel="Time", ylabel="Population", title=title)
    
    # Plot mean ± std
    plot!(p, stats.times, stats.P0_mean, 
          ribbon=stats.P0_std, fillalpha=0.3,
          label="ρ₀₀ (Ground)", linewidth=2, color=:blue)
    plot!(p, stats.times, stats.P1_mean, 
          ribbon=stats.P1_std, fillalpha=0.3,
          label="ρ₁₁ (Excited)", linewidth=2, color=:red)
    
    return p
end

"""
    plot_noise_trajectory(noise_sol; title="OU Noise Trajectory")

Plot a single noise trajectory.
"""
function plot_noise_trajectory(noise_sol; title="OU Noise Trajectory")
    p = plot(noise_sol.t, noise_sol[1,:], 
             label="ζ(t)", xlabel="Time", ylabel="Noise",
             title=title, linewidth=1.5)
    hline!([0], linestyle=:dash, color=:black, label="Mean", linewidth=1)
    return p
end

# ============================================================================
# Example Functions
# ============================================================================

"""
    example_pi_half()

Example: π/2 rotation starting from ground state (no noise).
"""
function example_pi_half()
    println("\n" * "="^60)
    println("Example: π/2 Rotation (X Gate)")
    println("="^60)
    
    # Simulate
    sol = simulate(ground_state(), 1.0, π/2)
    
    # Print results
    println("\nFinal state:")
    println("  ρ₀₀ = $(sol[1,end]) (expected ≈ 0.5)")
    println("  ρ₁₁ = $(sol[7,end]) (expected ≈ 0.5)")
    println("  Trace = $(sol[1,end] + sol[7,end]) (expected = 1.0)")
    
    # Plot
    p1 = plot_populations(sol, title="Population Dynamics")
    p2 = plot_coherences(sol, title="Coherence Evolution")
    display(plot(p1, p2, layout=(2,1), size=(800, 600)))
    
    return sol
end

"""
    example_pi()

Example: π rotation (complete population inversion).
"""
function example_pi()
    println("\n" * "="^60)
    println("Example: π Rotation (Population Inversion)")
    println("="^60)
    
    sol = simulate(ground_state(), 1.0, π)
    
    println("\nFinal state:")
    println("  ρ₀₀ = $(sol[1,end]) (expected ≈ 0.0)")
    println("  ρ₁₁ = $(sol[7,end]) (expected ≈ 1.0)")
    
    display(plot_populations(sol, title="π Rotation"))
    
    return sol
end

"""
    example_rabi()

Example: Multiple Rabi oscillations.
"""
function example_rabi()
    println("\n" * "="^60)
    println("Example: Rabi Oscillations (2π rotation)")
    println("="^60)
    
    sol = simulate(ground_state(), 2.0, 2π)
    
    display(plot_populations(sol, title="Rabi Oscillations"))
    
    return sol
end

"""
    example_no_noise()

Example: Clean π/2 rotation without noise.
"""
function example_no_noise()
    println("\n" * "="^60)
    println("Example: π/2 Rotation (No Noise)")
    println("="^60)
    
    sol = simulate(ground_state(), 1.0, π/2)
    
    println("\nFinal state:")
    println("  ρ₀₀ = $(sol[1,end]) (expected ≈ 0.5)")
    println("  ρ₁₁ = $(sol[7,end]) (expected ≈ 0.5)")
    
    display(plot_populations(sol, title="No Noise"))
    return sol
end

"""
    example_single_noise()

Example: π/2 rotation with single noise trajectory.
"""
function example_single_noise()
    println("\n" * "="^60)
    println("Example: π/2 Rotation (Single Noise Trajectory)")
    println("="^60)
    
    # Parameters
    tg = 1.0
    θ0 = π/2
    τ_c = 0.1  # Fast noise
    σ = 0.1    # 10% noise amplitude
    
    # Generate noise
    noise_sol = generate_ou_noise(τ_c, 0.0, σ, (0.0, tg), tg/1000)
    
    # Simulate with noise
    sol = simulate_with_noise(ground_state(), tg, θ0, noise_sol)
    
    println("\nFinal state:")
    println("  ρ₀₀ = $(sol[1,end])")
    println("  ρ₁₁ = $(sol[7,end])")
    
    # Plot
    p1 = plot_noise_trajectory(noise_sol)
    p2 = plot_populations(sol, title="With Noise")
    display(plot(p1, p2, layout=(2,1), size=(800, 600)))
    
    return sol, noise_sol
end

"""
    example_noise_ensemble()

Example: π/2 rotation with noise ensemble.
"""
function example_noise_ensemble()
    println("\n" * "="^60)
    println("Example: π/2 Rotation (Noise Ensemble)")
    println("="^60)
    
    # Parameters
    tg = 1.0
    θ0 = π/2
    τ_c = 0.1
    σ = 0.1
    n_samples = 50
    
    println("\nSimulating $n_samples trajectories...")
    solutions = simulate_noisy_ensemble(ground_state(), tg, θ0, τ_c, σ, n_samples)
    
    # Statistics
    stats = ensemble_statistics(solutions)
    println("\nFinal state statistics:")
    println("  ⟨ρ₀₀⟩ = $(stats.P0_mean[end]) ± $(stats.P0_std[end])")
    println("  ⟨ρ₁₁⟩ = $(stats.P1_mean[end]) ± $(stats.P1_std[end])")
    
    # Plot
    display(plot_ensemble_with_stats(solutions))
    
    return solutions
end

"""
    compare_noise_levels()

Compare effects of different noise amplitudes.
"""
function compare_noise_levels()
    println("\n" * "="^60)
    println("Example: Comparing Noise Levels")
    println("="^60)
    
    tg = 1.0
    θ0 = π/2
    τ_c = 0.1
    σ_values = [0.0, 0.05, 0.1, 0.2]
    n_samples = 30
    
    plots = []
    
    for σ in σ_values
        if σ == 0.0
            sol = simulate(ground_state(), tg, θ0)
            p = plot_populations(sol, title="σ = $σ (no noise)")
        else
            solutions = simulate_noisy_ensemble(ground_state(), tg, θ0, τ_c, σ, n_samples)
            p = plot_ensemble_with_stats(solutions, title="σ = $σ")
        end
        push!(plots, p)
    end
    
    display(plot(plots..., layout=(2,2), size=(1200, 900)))
end

# ============================================================================
# Print Help Message
# ============================================================================

println("""
╔════════════════════════════════════════════════════════════╗
║  Density Matrix + Noise Simulation Module Loaded!         ║
╚════════════════════════════════════════════════════════════╝

Clean Simulation (No Noise):
    sol = simulate(ground_state(), 1.0, π/2)
    plot_populations(sol)

Noisy Simulation (Single Trajectory):
    noise = generate_ou_noise(0.1, 0.0, 0.1, (0.0, 1.0), 0.001)
    sol = simulate_with_noise(ground_state(), 1.0, π/2, noise)

Noisy Ensemble:
    sols = simulate_noisy_ensemble(ground_state(), 1.0, π/2, 0.1, 0.1, 50)
    plot_ensemble_with_stats(sols)

Examples (No Noise):
    example_pi_half()            - π/2 rotation demo
    example_pi()                 - π rotation demo
    example_rabi()               - Rabi oscillations
    example_no_noise()           - Clean simulation

Examples (With Noise):
    example_single_noise()       - One noisy trajectory
    example_noise_ensemble()     - Ensemble average
    compare_noise_levels()       - Compare σ values
    
Key Functions:
    generate_ou_noise()          - Generate noise
    simulate_with_noise()        - Single noisy run
    simulate_noisy_ensemble()    - Multiple noisy runs
    ensemble_statistics()        - Get mean/std
    plot_ensemble_with_stats()   - Plot ensemble
""")