using ProgressLogging ,DifferentialEquations, Plots, Statistics, DifferentialEquations.EnsembleAnalysis, DiffEqNoiseProcess 

function OrsteinUhlebeckAddaptive(θ, μ, σ, u0, tspan, saveint, solver, W)
    
    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    
    prob = SDEProblem(f, g, u0, tspan, noise = W)
    sol = solve(prob, solver, saveat = saveint)
    
    return sol.u

end

function AverageValueAddaptive(θ2, μ2, σ2, u02, tspan2, saveint2, solver2, samp)

    f(u, p, t) = θ2 * (μ2 - u)
    g(u, p, t) = σ2
    prob = SDEProblem(f, g, u02, tspan2)
    ensembleprob = EnsembleProblem(prob)
    sim = solve(ensembleprob, saveat=saveint2, solver2,EnsembleThreads(), trajectories = samp)
    sol = sim.u
    covs = cov(sol)

    return mean(sol), covs

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
    x = Int(samp3/100)

    for i in x:x:samp3

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

function OrsteinUhlebeckFixed(θ, μ, σ, u0, tspan, dt1, solver, reltol1, abstol1)
    
    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, dt = dt1, abstol = abstol1, reltol = reltol1)
    
    return sol.u

end

function AverageValueFixed(θ2, μ2, σ2, u02, tspan2, dt2, solver2, samp, reltol2, abstol2)

    f(u, p, t) = θ2 * (μ2 - u)
    g(u, p, t) = σ2
    
    prob = SDEProblem(f, g, u02, tspan2)
    ensembleprob = EnsembleProblem(prob)
    sim = solve(ensembleprob, dt=dt2, solver2,EnsembleThreads(), trajectories = samp, abstol=abstol2, reltol=reltol2)
    sol = sim.u
    covs = cov(sol)
    return mean(sol), covs
end

function sampletimeFixed(θ3, μ3, σ3, u03, tspan3, dt3, solver3, samp2, reltol3, abstol3)
    
    elapsed_time = @elapsed begin
    
        sol = AverageValueFixed(θ3, μ3, σ3, u03, tspan3, dt3, solver3, samp2, reltol3, abstol3)
    
        t = range(0, stop = tspan3[2], step = dt3)
    
    end

    s = elapsed_time

    return s

end

function timeperSampleFixed(θ4, μ4, σ4, u04, tspan4, dt4, solver4, samp3, reltol4, abstol4)

    tp = []
    x = Int(samp3/100)

    for i in x:x:samp3

        ts = sampletimeFixed(θ4, μ4, σ4, u04, tspan4, dt4, solver4, i, reltol4, abstol4)

        push!(tp,ts)

    end

    return tp

end

function timeperSampleAverageFixed(θ5, μ5, σ5, u05, tspan5, dt5, solver5, samp4, amount, reltol5, abstol5)

    alltps = []

    for i in 1:amount

        ts = timeperSampleFixed(θ5, μ5, σ5, u05, tspan5, dt5, solver5, samp4, reltol5, abstol5)

        push!(alltps,ts)

    end

    return sum(alltps, dims=1)/amount

end

function OrnsteinUhlenbeckAvg(θ, μ, u0, t)
    
    u0*exp(-θ*t) + μ*(1-exp(-θ*t))
     
end

function OrnsteinUhlenbeckAvgValues(θ, μ, u0, t)

    values = []

    for i in t

        x = OrnsteinUhlenbeckAvg(θ, μ, u0, i)

        push!(values, x)
    
    end

    return values

end

function ErrorCalculator(x,y)

    errors = []

    for i in 1:lastindex(x)

        z = abs((y[i]-x[i])/x[i])*100

        push!(errors,z)
    end

    return errors

end

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

soln, cov1 = AverageValueAddaptive(θ, μ, σ, u0, tspan, saveint, LambaEM(), samples)
soln2, cov2 = AverageValueFixed(θ, μ, σ, u0, tspan, dt2, EM(), samples, reltol, absol)
soln3, cov3 = AverageValueAddaptive(θ, μ, σ, u0, tspan, saveint, LambaEulerHeun(), samples)
soln4, cov4 = AverageValueFixed(θ, μ, σ, u0, tspan, dt2, EulerHeun(), samples, reltol, absol)


tims =  timeperSampleAverageAddaptive(θ, μ, σ, u0, tspan, saveint, LambaEM(), samples, amount)
tims2 = timeperSampleAverageFixed(θ, μ, σ, u0, tspan, dt2, EM(), samples, amount, reltol, absol)
tims3 =  timeperSampleAverageAddaptive(θ, μ, σ, u0, tspan, saveint, LambaEulerHeun(), samples, amount)
tims4 = timeperSampleAverageFixed(θ, μ, σ, u0, tspan, dt2, EulerHeun(), samples, amount, reltol, absol)

t = 10.0
dt1 = saveint
arr = range(0, 10, 101)
arr1 = range(0, 10, Int(tf/dt2 + 1))

soln5 = OrnsteinUhlenbeckAvgValues(θ, μ, u0, arr)
soln6 = OrnsteinUhlenbeckAvgValues(θ, μ, u0, arr1)

error1 = ErrorCalculator(soln5, soln)
error2 = ErrorCalculator(soln6, soln2)
error3 = ErrorCalculator(soln5, soln3)
error4 = ErrorCalculator(soln6, soln4)

arr2 = 1:lastindex(tims)

y = Int(samples/100)

x = y:y:samples

p3 = plot(x,tims, xlabel = "samples", ylabel = "Time", label = "ito addaptive step size") 
plot!(x,tims2, xlabel = "samples", ylabel = "Time", label = "ito fixed step size") 
plot!(x,tims3, xlabel = "samples", ylabel = "Time", label = "Euler-Heun addaptive step size") 
plot!(x,tims4, xlabel = "samples", ylabel = "Time", label = "Euler-Heun fixed step size")

p1 = plot(soln, xlabel = "Time", ylabel = "Y values", label = "ito addaptive step size") 
plot!(soln2, xlabel = "Time", ylabel = "Y values", label = "ito fixed step size") 
plot!(soln3, xlabel = "Time", ylabel = "Y values", label = "Euler-Heun addaptive step size") 
plot!(soln4, xlabel = "Time", ylabel = "Y values", label = "Euler-Heun fixed step size") 
plot!(soln5, xlabel = "Time", ylabel = "Y values", label = "OrnsteinUhlenbeckAvgValues")

p2 = plot(arr,error1, xlabel = "Time", ylabel = "Error", label = "ito addaptive step size") 
plot!(arr1,error2, xlabel = "Time", ylabel = "Error", label = "ito fixed step size") 
plot!(arr,error3, xlabel = "Time", ylabel = "Error", label = "Euler-Heun addaptive step size") 
plot!(arr1,error4, xlabel = "Time", ylabel = "Error", label = "Euler-Heun fixed step size")

p4 = plot(heatmap(cov1))
plot!(heatmap(cov2))
plot!(heatmap(cov3))
plot!(heatmap(cov4))

display(p1)
display(p2)
display(p3)
display(p4)