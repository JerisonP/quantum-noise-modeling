using DiffEqNoiseProcess, ProgressLogging ,DifferentialEquations, Plots, Statistics, SignalAnalysis

#construct noise 

θ = 1.0
μ = 0
σ = sqrt(2)

f(u, p, t) = θ * (μ - u)
g(u, p, t) = σ
W = OrnsteinUhlenbeckProcess(1.0 , 0, sqrt(2) , 0.0, 1)

prob = SDEProblem(f, g, u0, tspan, noise = W)
sol = solve(prob, LambaEulerHeun(), saveat = 0.1)

plot(sol)