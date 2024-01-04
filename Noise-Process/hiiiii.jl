using DifferentialEquations
using Plots
using Statistics
using SignalAnalysis

# Define the parameters
sampling_rate = 1000  # Sampling rate (Hz)
total_time = 10        # Total time duration (seconds)
noise_variance = 0.5   # Variance of the noise

# Generate time array
tspan = (0.0, total_time)
t = range(tspan[1], tspan[2], length = total_time * sampling_rate)

# Define the drift and noise functions for the process
function linear_drift!(du, u, p, t)
    du = 1.0  # Change this function as per your model
end

function noise_func!(du, u, p, t)
    du .= p[1] * white_noise(u, p, t)  # Adding colored noise
end

# Set up the colored noise process
p = [noise_variance]
prob = SDEProblem(linear_drift!, noise_func!, ones(1), tspan, p)

# Solve the SDE to generate the noise
sol = solve(prob, SRIW1(), noise = true, saveat = t)

# Extract the generated noise from the solution
noise = sol[end]

# Plot the generated colored noise
plot(t, noise, xlabel = "Time", ylabel = "Amplitude", title = "Colored Noise")
