using Plots
using Statistics
using StatsPlots

using Printf # for formating floats as strings
using DSP #for conv and bartlett smoothing

include("./coorddecent.jl")

# Import data
include("./dataimport.jl")
filename = "Lee et al 2021 distributions.xlsx"
distributions = read_distribution_data(filename);
D = to_tensor(distributions)

# Preprocessing helpers to smooth out distribution
smooth(v,n::Integer;rescale=true) = conv(v,bartlett(2*n+1))[n+1:end-n] / (rescale ? n : 1)
coarsen(v, n) = smooth(v,n,rescale=false)[1:n÷2:end]/2
function normalize_fiber_sum!(D, dims)
    ds = eachslice(D, dims=dims)
    ds ./= (sum.(ds))
end

# Ensure every distribution has area 1
normalize_fiber_sum!(D, (2,3));

# Look at distributions
k = 2
n = 250
(; name, scale) = distributions[k]
name=short_name(name)

plot(scale, [d for d ∈ eachcol(D[:,:,k])],
    title="Raw $name Distribution for Each Rock",
    yaxis="mass",
    xaxis= (name=="Age") ? "Ma" : "ppm")

plot(scale, [smooth(d, n) for d ∈ eachcol(D[:,:,k])],
    title="$name Distribution for Each Rock smoothed to n=$n",
    yaxis="mass",
    xaxis= (name=="Age") ? "Ma" : "ppm")

plot(scale[1:n÷2:end], [coarsen(d, n) for d ∈ eachcol(D[:,:,k])],
    title="$name Distribution for Each Rock coarsened to n=$n",
    yaxis="mass",
    xaxis= (name=="Age") ? "Ma" : "ppm")

# Smooth out distributions according to:
# length(density vector) = n * √(# of grain samples per rock)
# 2000 ≈ 250 * √(60)
# This gives us histogram-like bars where the with is around √ the number of samples
downsample = true
pre_process = downsample ? coarsen : smooth
D_processed = mapslices(d -> pre_process(d, n), D, dims=1);

# Perform Decomposition
T = permutedims(D_processed,(2,3,1));
normalize_fiber_sum!(T, (1,2));
fiber_sums_normalized = all(sum.(eachslice(T,dims=(1,2))) .≈ 1.0)
#findall(0.99 .> sum.(eachslice(T,dims=(1,2))) )
number_of_factors = 2
C, F, error, norm_grad = coorddecent(T, number_of_factors, maxiter=10000, tol=2e-5); #1e-6 for first data set with error convergence

# Rescale to interpret results
fiber_sums = sum.(eachslice(F,dims=(1,2)))
avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
F = avg_factor_sums^(-1) * F # TODO make more accurate scaling
C = C * avg_factor_sums

# Plots the resulting coefficients
heatmap(C,
    yticks=(eachindex(C[:,1]), ["rock $i" for i ∈ eachindex(C[:,1])]),
    xticks=eachindex(C[1,:]),
    xlabel="Factor",
    title="Amount of Factors in Rock Samples"
)

# For each element, compare distributions between factors
#element_list = short_name.(distributions.names)
for (j, d) ∈ enumerate(distributions)
    F_slice = F[:,j,:]
    (; name, scale) = d
    name = short_name(name)
    xlabel = (name == "Age") ?  "Ma" : "ppm"
    xscale = downsample ? scale[1:n÷2:end] : scale
    p = heatmap(xscale, (eachindex(F_slice[:,1])), F_slice,
        xlabel = xlabel,
        ylabel = "Factor",
        title = "Distribution of $name in Factors",
        yticks = eachindex(F_slice[1,:])
    )
    display(p)
end

element_list = short_name.(distributions.names)
# Show distribution of different elements within a factor
for (i, F_slice) ∈ enumerate(eachslice(F,dims=1))
    p = heatmap(F_slice,
        yticks=(eachindex(F_slice[:,1]), element_list),
        xticks=([1, length(F_slice[1,:])],["0", "1"]),
        xlabel="Normalized Range of Values",
        title = "Distribution of Elements in Factor $i"
        )
    display(p)
end