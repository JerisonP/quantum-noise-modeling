using Pkg
Pkg.add("Distributions")

using Distributions
using Plots

# Create a random sigma values to plot
sigmas = []

for i in 1000:2000
    append!(sigmas, sqrt(2)/i)

end

random = rand(1:1000)

println(sigmas[random])

d = Normal(0, sigmas[random])
    
# Draw 100 samples from the distribution
x = rand(d, 10000)
    
# Display the plots
plot(histogram(x, title="sigma = $sigmas[random]", label=""))
