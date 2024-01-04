include("C:\\Users\\whitl\\OneDrive\\Desktop\\My Stuff\\Programing\\Research\\OUPLibrary.jl")

θ = 1.0
μ = 0
σ = sqrt(2)

u0 = 10.0  
tf = 10

tspan = (0.0, tf)  # Time interval from 0.0 to 10.0

samples = 100000
saveint = 0.5

solver = LambaEM()
solver2 = LambaEulerHeun()
solver3 = EM()

soln = AverageValueAddaptive(θ, μ, σ, u0, tspan, saveint, solver, samples)
soln2 = AverageValueAddaptive(θ, μ, σ, u0, tspan, saveint, solver2, samples)
soln3 = AverageValueFixed(θ, μ, σ, u0, tspan, dt, reltol, abstol, solver3, samples)
t = TimeArry(tspan, saveint)

plot(t,soln, xlabel = "Time", ylabel = "u", label = "ito addaptive step size")
plot!(arr1,soln3, xlabel = "Time", ylabel = "u", label = "ito fixed step size")
plot!(arr,soln2, xlabel = "Time", ylabel = "u", label = "modified Euler-Heun method adaptive step size")
plot!(arr,OrnsteinUhlenbeckAvg.(θ, μ, u0, arr), xlabel = "Time", ylabel = "u", title = "Ornstein-Uhlenbeck Process")