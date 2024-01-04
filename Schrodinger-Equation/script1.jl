using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")

using DifferentialEquations
using Plots

# Define the ODEs as a function
function f(du,u,p,t)
    O = pi/2
    tf = 1000
    ft = O/tf*(1-cos(2*pi*t/tf))
    # Separate real and imaginary parts
    Ar, Ai, Br, Bi = u
    # Define the derivatives for the real and imaginary parts
    du[1] = -0.5 * ft * Bi
    du[2] = 0.5 * ft * Br
    du[3] = -0.5 * ft * Ai
    du[4] = 0.5 * ft * Ar
end

tfinal = 1000 # Set the final time you want to solve the ODEs up to
u0 = [1.0, 0.0, 0.0, 0.0] # Set the initial conditions for the real and imaginary parts
tspan = (0.0, tfinal) # Set the time span to solve the ODEs over

# Create an ODEProblem object
prob = ODEProblem(f, u0, tspan)

# Solve the ODEs using the Vern9 solver
sol = solve(prob, Feagin14(), reltol=1e-18, abstol=1e-16)

# Combine the real and imaginary parts of the solution
A = sol[1,:] .+ im .* sol[2,:]
B = sol[3,:] .+ im .* sol[4,:]

# Plot the solution
plot(sol.t, [abs2.(A) abs2.(B)], label=["|A|^2" "|B|^2"], xlabel="Time", ylabel="Magnitude squared")
# Save the plot to a file
#savefig("plot.png")
