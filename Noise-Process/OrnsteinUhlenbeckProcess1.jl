using DifferentialEquations, Plots

function OrnsteinUhlenbeck(θ, μ, σ, u0, tspan, saveint)

    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ

    prob = SDEProblem(f, g, u0, tspan)
    #sol = solve(prob, SRIW1(), saveat = saveint)
    sol = solve(prob, LambaEulerHeun(), saveat = saveint)

    return sol.u

end

function AverageValue(θ, μ, σ, u0, tspan, samples, saveint)

    ODEValTot = []

    for j in 1:samples

        ODEvalues = OrnsteinUhlenbeck(θ, μ, σ, u0, tspan, saveint)

        push!(ODEValTot,ODEvalues)

    end

    return sum(ODEValTot, dims=1)/samples

end

function OrnsteinUhlenbeckAvg(θ, μ, u0, t)
    
   u0*exp(-θ*t) + μ*(1-exp(-θ*t))
    
end
