import Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
using DifferentialEquations, Plots

function OrsteinUhlebeckAddaptive(θ, μ, σ, u0, tspan, saveint, solver)
    
    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, saveat = saveint)
    
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

function timeperSampleAddaptive(θ4, μ4, σ4, u04, tspan4, saveint4, solver4, samp3)

    tp = []

    for i in 1:samp3

        ts = sampletimeAddaptive(θ4, μ4, σ4, u04, tspan4, saveint4, solver4, i)

        push!(tp,ts)

    end

    return tp

end

function timeperSampleAverageAddaptive(θ5, μ5, σ5, u05, tspan5, saveint5, solver5, samp4, amount)

    alltps = []

    for i in 1:amount

        ts = timeperSampleAddaptive(θ5, μ5, σ5, u05, tspan5, saveint5, solver5, samp4)

        push!(alltps,ts)

    end

    return sum(alltps, dims=1)/amount

end

function OrsteinUhlebeckFixed(θ, μ, σ, u0, tspan, dt1, solver)
    
    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, dt = dt1)
    
    return sol.u

end

function AverageValueFixed(θ2, μ2, σ2, u02, tspan2, dt2, solver2, samp)

    ODEValTot = []

    for j in 1:samp

        ODEvalues = OrsteinUhlebeckFixed(θ2, μ2, σ2, u02, tspan2, dt2, solver2)
    
        push!(ODEValTot,ODEvalues)
    
    end

    Val = sum(ODEValTot, dims=1)/samp

    return Val[1]

end

function sampletimeFixed(θ3, μ3, σ3, u03, tspan3, dt3, solver3, samp2)
    
    elapsed_time = @elapsed begin
    
        sol = AverageValueFixed(θ3, μ3, σ3, u03, tspan3, dt3, reltol3, abstol3, solver3, samp2)
    
        t = range(0, stop = tspan3[2], step = dt3)
    
    end

    s = elapsed_time

    return s

end

function timeperSampleFixed(θ4, μ4, σ4, u04, tspan4, dt4, solver4, samp3)

    tp = []

    for i in 1:samp3

        ts = sampletimeFixed(θ4, μ4, σ4, u04, tspan4, dt4, solver4, i)

        push!(tp,ts)

    end

    return tp

end

function timeperSampleAverageFixed(θ5, μ5, σ5, u05, tspan5, dt5, solver5, samp4, amount)

    alltps = []

    for i in 1:amount

        ts = timeperSampleFixed(θ5, μ5, σ5, u05, tspan5, dt5, solver5, samp4)

        push!(alltps,ts)

    end

    return sum(alltps, dims=1)/amount

end

function OrnsteinUhlenbeckAvg(θ, μ, u0, t)
    
    u0*exp(-θ*t) + μ*(1-exp(-θ*t))
     
end