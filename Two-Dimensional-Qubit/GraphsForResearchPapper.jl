using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Interact")
Pkg.add("Graphs")
Pkg.add("Plots")

using DifferentialEquations
using Interact
using Graphs
using Plots

default(fontfamily="sans-serif")

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

t = (1:1000)

# Combine the real and imaginary parts of the solution
A = sol(t, idxs = 1) .+ im .* sol(t, idxs = 2)
B = sol(t, idxs = 3) .+ im .* sol(t, idxs = 4)

C = abs2.(A)
D = abs2.(B)
tScaled = t/tfinal

file = open("output.txt", "w")
write(file, "tScaled  C  D\n")

for i in 1:length(C)
    a = tScaled[i]
    b = C[i]
    c = D[i]
    write(file, "$a  $b  $c\n")
end

close(file)
