using ProgressLogging ,DifferentialEquations, Plots, Statistics, DifferentialEquations.EnsembleAnalysis, DiffEqNoiseProcess 

θ = 1.0
μ = 0
σ = sqrt(2)

u0 = 10.0  
tf = 10

tspan = (0.0, tf)  # Time interval from 0.0 to 10.0

samples = 1000
saveint = 0.1
dt2 = 0.1

amount = 2

reltol = absol = 1.0*10^(-6)

function OrsteinUhlebeckAddaptive(θ, μ, σ, u0, tspan, saveint, solver, W)

    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    W = OrnsteinUhlenbeckProcess(θ,μ,σ,0.0,1.0)
    prob = SDEProblem(f, g, u0, tspan, noise = W)
    solve(prob, LambaEM(), saveat = saveint)
    return sol.u
end

function AverageValueAddaptive(θ2, μ2, σ2, u02, tspan2, saveint2, solver2, samp)

    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    W = OrnsteinUhlenbeckProcess(3.0,0.0,2.7*σ,0.0,0.0)
    prob = SDEProblem(f, g, u0, tspan, noise = W)
    ensembleprob = EnsembleProblem(prob)
    sim = solve(ensembleprob, saveat=saveint2, solver2,EnsembleThreads(), trajectories = samp)
    sol = sim.u
    covs = cov(sol)

    return mean(sol), covs

end

function AverageValueAddaptive1(θ2, μ2, σ2, u02, tspan2, saveint2, solver2, samp)

    f(u, p, t) = θ2 * (μ2 - u)
    g(u, p, t) = σ2
    prob = SDEProblem(f, g, u02, tspan2)
    ensembleprob = EnsembleProblem(prob)
    sim = solve(ensembleprob, saveat=saveint2, solver2,EnsembleThreads(), trajectories = samp)
    sol = sim.u
    covs = cov(sol)

    return mean(sol), covs

end

x, y = AverageValueAddaptive(θ, μ, σ, u0, tspan, saveint, LambaEM(), samples)

x_1, y_1 = AverageValueAddaptive1(θ, μ, σ, u0, tspan, saveint, LambaEM(), samples)

plot(x)
plot!(x_1)
