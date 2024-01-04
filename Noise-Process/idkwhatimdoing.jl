using ProgressLogging ,DifferentialEquations, Plots, Statistics, SignalAnalysis

function OrsteinUhlebeckAddaptive(θ, μ, σ, u0, tspan, saveint, solver)
    
    f(u, p, t) = u
    g(v, p, t) = θ * (μ - v)
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, saveat = saveint, progress = true)
    
    return sol.u

end

function AverageValueAddaptive(θ2, μ2, σ2, u02, tspan2, saveint2, solver2, samp)

    ODEValTot = []

    for j in 1:samp

        ODEvalues = OrsteinUhlebeckAddaptive(θ2, μ2, σ2, u02, tspan2, saveint2, solver2)
    
        push!(ODEValTot,ODEvalues)
    
    end

    Val = sum(ODEValTot, dims=1)/samp

    return Val[1]

end

function TimeArry(tspan, dt)
    
    return range(tspan[1], tspan[2], Int((tspan[2]-tspan[1])/dt + 1))

end

function sampletimeAddaptive(θ3, μ3, σ3, u03, tspan3, saveint3, solver3, samp2)
    
    elapsed_time = @elapsed begin
    
        sol = AverageValueAddaptive(θ3, μ3, σ3, u03, tspan3, saveint3, solver3, samp2)
    
        t = range(0, stop = tspan3[2], step = saveint3)
    
    end

    s = elapsed_time

    return s

end

θ = 1.0
μ = 0
σ = sqrt(2)

u0 = 10.0  
tf = 10

tspan = (0.0, tf)  # Time interval from 0.0 to 10.0

samples = 10
saveint = 0.1
dt2 = 0.1

soln3 = AverageValueAddaptive(θ, μ, σ, u0, tspan, saveint, LambaEulerHeun(), samples)

t = 10.0
dt1 = saveint
arr = range(0, tf, Int(tf/dt1 + 1))

plot(arr,soln3, xlabel = "Time", ylabel = "Y values", label = "Euler-Heun addaptive step size")
