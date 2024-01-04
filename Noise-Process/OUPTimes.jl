import Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
Pkg.add("TickTock")
using DifferentialEquations, Plots

function OrsteinUhlebeck(θ, μ, σ, u0, tspan, saveint, solver)
    
    f(u, p, t) = θ * (μ - u)
    g(u, p, t) = σ
    
    prob = SDEProblem(f, g, u0, tspan)
    sol = solve(prob, solver, saveat = saveint)
    
    return sol.u

end

function AverageValue(θ2, μ2, σ2, u02, tspan2, saveint2, solver2, samp)

    ODEValTot = []

    for j in 1:samp

        ODEvalues = OrsteinUhlebeck(θ2, μ2, σ2, u02, tspan2, saveint2, solver2)
    
        push!(ODEValTot,ODEvalues)
    
    end

    return sum(ODEValTot, dims=1)/samp

end

function sampletime(θ3, μ3, σ3, u03, tspan3, saveint3, solver3, samp2)
    
    elapsed_time = @elapsed begin
    
        sol = AverageValue(θ3, μ3, σ3, u03, tspan3, saveint3, solver3, samp2)
    
        t = range(0, stop = tspan3[2], step = saveint3)
    
    end

    s = elapsed_time

    return s

end

function timeperSample(θ4, μ4, σ4, u04, tspan4, saveint4, solver4, samp3)

    tp = []

    for i in 1:samp3

        ts = sampletime(θ4, μ4, σ4, u04, tspan4, saveint4, solver4, i)

        push!(tp,ts)

    end

    popfirst!(tp)

    return tp

end


#=
function timeperSampleAverage(θ5, μ5, σ5, u05, tspan5, saveint5, solver5, samp4, amount)

    alltps = []

    for i in 1:amount

        ts = timeperSample(θ5, μ5, σ5, u05, tspan5, saveint5, solver5, samp4)

        push!(alltps,ts)

    end

    return sum(alltps, dims=1)/amount

end
=#

θ = 1.0
μ = 0
σ = sqrt(2)

u0 = 10.0  
tf = 10

tspan = (0.0, tf)  # Time interval from 0.0 to 10.0

samp = 100

saveint = 0.25

solver2 = LambaEM()
solver4 = LambaEulerHeun()
solver5 = RKMil()
solver6 = RKMilCommute()
solver7 = RKMilGeneral()
solver14 = SRA()
solver15 = SRI()
solver16 = SRIW1()
solver17 = SRIW2()
solver18 = SOSRI()
solver19 = SOSRI2()
solver20 = SRA1()
solver21 = SRA2()
solver22 = SRA3()
solver23 = SOSRA()
solver24 = SOSRA2()
solver30 = ImplicitEM()
solver31 = STrapezoid()
solver32 = SImplicitMidpoint()
solver33 = ImplicitEulerHeun()
solver34 = ImplicitRKMil()
solver35 = ISSEM()
solver36 = ISSEulerHeun()
solver39 = DRI1()
solver41 = RI1()
solver42 = RI3()
solver43 = RI5()
solver44 = RI6()
solver46 = RDI2WM()
solver47 = RDI3WM()
solver48 = RDI4WM()
solver59 = DRI1()

x2 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver2, samp)
x4 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver4, samp)
x5 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver5, samp)
x6 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver6, samp)
x7 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver7, samp)
x14 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver14, samp)
x15 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver15, samp)
x16 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver16, samp)
x17 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver17, samp)
x18 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver18, samp)
x19 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver19, samp)
x20 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver20, samp)
x21 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver21, samp)
x22 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver22, samp)
x23 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver23, samp)
x24 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver24, samp)
x30 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver30, samp)
x31 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver31, samp)
x32 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver32, samp)
x33 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver33, samp)
x34 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver34, samp)
x35 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver35, samp)
x36 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver36, samp)
x39 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver39, samp)
x41 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver41, samp)
x42 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver42, samp)
x43 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver43, samp)
x44 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver44, samp)
x46 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver46, samp)
x47 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver47, samp)
x48 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver48, samp)
x59 = timeperSample(θ, μ, σ, u0, tspan, saveint, solver59, samp)

plot( x16, xlabel = "samples", ylabel = "time16", title = "Ornstein-Uhlenbeck Process")
