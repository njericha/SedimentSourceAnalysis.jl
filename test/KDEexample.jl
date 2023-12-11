###############
# This file compares the histogram to a KDE with box, Gaussian, and Dirac kernels
###############

using Random
using Plots
using KernelDensity

#using Pkg
#Pkg.add("Distributions")
using Distributions

# Make Distribution
"""
Returns a sample Z from a distribution D where
- Z ∼ N(2, 1) with probability 0.2
- Z ∼ N(-4, 2) with probability 0.8
"""
function D(p=0.2)
    flip = rand()
    if flip <= p
        return randn() + 2 # Z ∼ N(2, 1)
    else
        return 2*randn() - 4 # Z ∼ N(-4, 2)
    end
end

MixtureModel(Normal[
   Normal(2, 1),
   Normal(-4, 2)], [0.2, 0.8])

# Sample distribution N times
N = 100
samples = [D() for _ in 1:N]

# Plot the histogram and samples
bins = ceil(Int,sqrt(N))
histogram(samples, bins=bins, normalize=true; label = "histogram") |> display
scatter!(samples, zero(samples); markershape=:vline, markersize=20, label = "Samples", color=:black)

# Plot the Gaussian KDE
f = kde(samples)
plot!(f.x, f.density; linewidth=7, color=:orange, label = "Gaussian KDE")

# Plot the Uniform (box car) KDE
g = kde(samples, kernel=Uniform)
plot!(g.x, g.density; linewidth=4, label = "Uniform KDE", color=:purple, alpha=0.8)

# Make density function
gaussian(μ=0, σ=1) = x -> pdf(Normal(μ, σ), x)

# Plot the true density
true_density(x) = 0.2*(gaussian(2, 1))(x) + 0.8*(gaussian(-4, 2))(x)
plot!(f.x, true_density.(f.x);
    linewidth=3,
    linestyle=:dash,
    label = "true density",
    xlabel = "value",
    ylabel = "density",
    color=:black,
    alpha=0.6,
    ylims = [0, 0.2])
