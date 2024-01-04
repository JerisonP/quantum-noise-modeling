using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Interpolations")
Pkg.add("Distributions")
Pkg.add("Plots")

using DifferentialEquations
using Interpolations
using Distributions
using Plots

function omega(σ1, range1 = 10000)

    # Create distribution
    d = Normal(0, σ1) 

    # Draw 10000 samples from the distribution
    return rand(d, range1)

end

function ODESolver(w = 0, tf1 = 1000, u0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan = (0.0, tf), O = pi/2)

    function f(du,u,p,t)

        ft = O/tf1*(1-cos(2*pi*t/tf1))

        # Separate real and imaginary parts
        Ar, Ai, Br, Bi, Cr, Ci, Dr, Di = u
        
        # Define the derivatives for the real and imaginary parts
        du[1] = -0.5 * ft * (Ci-Bi) # Ar
        du[2] = 0.5 * ft * (Cr-Br) # Ai
        du[3] = -w * Bi - 0.5 * ft * (Di-Ai) # Br
        du[4] = w * Br + 0.5 * ft * (Dr-Ar) # Bi
        du[5] = w * Ci - 0.5 * ft * (Ai-Di) # Cr
        du[6] = -w * Cr + 0.5 * ft * (Ar-Dr) # Ci
        du[7] = -0.5 * ft * (Bi-Ci) # Dr
        du[8] = 0.5 * ft * (Br-Cr) # Di

    end

    # Create an ODEProblem object
    prob = ODEProblem(f,u0,tspan)

    # Solve the ODEs using the Vern9 solver
    sol = solve(prob, Feagin14(), reltol=1e-15, abstol=1e-15, saveat = 1)
    
    return sol

end

function AverageValue(x, tf2, range2)

    ODEValTot = []

    for j in 1:1:length(x)

        ODEvalues = ODESolver(x[j],tf2)
    
        push!(ODEValTot,ODEvalues)
    
    end

    return sum(ODEValTot, dims=1)/range2

end