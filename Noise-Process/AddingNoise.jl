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


function ValueExtractor(Averageρ)
    
    Arv = []
    Aiv = []
    Brv = []
    Biv = []
    Crv = []
    Civ = []
    Drv = []
    Div = []

    for i in 1:tf+1
        push!(Arv, Averageρ[1][i][1])
        push!(Aiv, Averageρ[1][i][2])
        push!(Brv, Averageρ[1][i][3])
        push!(Biv, Averageρ[1][i][4])
        push!(Crv, Averageρ[1][i][5])
        push!(Civ, Averageρ[1][i][6])
        push!(Drv, Averageρ[1][i][7])
        push!(Div, Averageρ[1][i][8])
    end

    Average = [Arv Aiv Brv Biv Crv Civ Drv Div]

    return Average'

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
    
tf = 1000
range = 10000
pro = sqrt(2)
σ = pro/tf

x = omega(σ, range)

Averageρ = AverageValue(x, tf, range )

ODETrueValues = ODESolver(0,tf)

t = collect(0:tf)/tf

ReshapedAverageρ = ValueExtractor(Averageρ)

Arv = ReshapedAverageρ[1,:]
Aiv = ReshapedAverageρ[2,:]
Brv = ReshapedAverageρ[3,:]
Biv = ReshapedAverageρ[4,:]
Crv = ReshapedAverageρ[5,:]
Civ = ReshapedAverageρ[6,:]
Drv = ReshapedAverageρ[7,:]
Div = ReshapedAverageρ[8,:]

Ar = ODETrueValues[1,:]
Ai = ODETrueValues[2,:]
Br = ODETrueValues[3,:]
Bi = ODETrueValues[4,:]
Cr = ODETrueValues[5,:]
Ci = ODETrueValues[6,:]
Dr = ODETrueValues[7,:]
Di = ODETrueValues[8,:]



A = Ar .+ im .* Ai
B = Br .+ im .* Bi
C = Cr .+ im .* Ci
D = Dr .+ im .* Di

Av = Arv .+ im .* Aiv
Bv = Brv .+ im .* Biv
Cv = Crv .+ im .* Civ
Dv = Drv .+ im .* Div

print("σ*tf = ")
print(σ*tf)
plot(t, [abs.(Av) abs.(Dv)], label=["P_10" "P_01"], xlabel="Time t/tf", linecolor=[:red :blue], ylabel="Probability")
plot!(t, [abs.(A) abs.(D)], label=["True P_10" "True P_01"], linecolor=[:red :blue], linestyle=:dash)