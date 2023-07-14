#=
This script uses the distributions made by Joel
=#

using Plots
using Statistics
using StatsPlots

using Printf # for formating floats as strings
using DSP #for conv and bartlett smoothing

include("./coorddecent.jl")

# Import data
include("./dataimport.jl")
filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022_DISTRIBUTIONS.xlsx"
distributions = read_distribution_data(filename);
D = to_tensor(distributions);

# Preprocessing helpers to smooth out distribution
smooth(v,n::Integer;rescale=true) = conv(v,bartlett(2*n+1))[n+1:end-n] / (rescale ? n : 1)
coarsen(v, n) = smooth(v,n,rescale=false)[1:n÷2:end]/2
function normalize_fiber_sum!(D, dims)
    ds = eachslice(D, dims=dims)
    ds ./= (sum.(ds))
end

# Ensure every distribution has area 1
normalize_fiber_sum!(D, (2,3))

# Look at distributions
k = 3
n = 100
(; name, scale) = distributions[k]
#name=short_name(name)

sample_name = "Sink"
plot(scale, [d for d ∈ eachcol(D[:,:,k])],
    title="Raw $name Distribution for Each $sample_name",
    yaxis="mass",
    xaxis= (name=="Ages") ? "Ma" : "ppm")

plot(scale, [smooth(d, n) for d ∈ eachcol(D[:,:,k])],
    title="$name Distribution for Each $sample_name smoothed to n=$n",
    yaxis="mass",
    xaxis= (name=="Ages") ? "Ma" : "ppm")

plot(scale[1:n÷2:end], [coarsen(d, n) for d ∈ eachcol(D[:,:,k])],
    title="$name Distribution for Each $sample_name coarsened to n=$n",
    yaxis="mass",
    xaxis= (name=="Ages") ? "Ma" : "ppm")

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
number_of_factors = 3
C, F, error, norm_grad = coorddecent(T, number_of_factors, maxiter=10000, tol=1e-5); #1e-6 for first data set with error convergence

# Rescale to interpret results
fiber_sums = sum.(eachslice(F,dims=(1,2)))
avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
F = avg_factor_sums^(-1) * F # TODO make more accurate scaling
C = C * avg_factor_sums

# Plot gradient convergence
plot(norm_grad ./ norm_grad[1],
    yaxis=:log10,
    ylabel="2-norm of Full Gradient",
    title="Convergence of Gradient",
    xlabel="iteration #")

# Plots the resulting coefficients
heatmap(C,
    yticks=(eachindex(C[:,1]), ["$sample_name $i" for i ∈ eachindex(C[:,1])]),
    xticks=eachindex(C[1,:]),
    xlabel="Factor",
    title="Amount of Factors in $sample_name Samples"
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

element_list = distributions.names
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


#= True Source proportions =#
source_filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
source_proportions = XLSX.readdata(source_filename,"Source proportions","B2:D21")
source_proportions ./= 75 # 75 Grains in each sink

C_true = @view source_proportions[:,[1,2,3]]

heatmap(C_true,
    yticks=(eachindex(C_true[:,1]), ["$sample_name $i" for i ∈ eachindex(C_true[:,1])]),
    xticks=eachindex(C_true[1,:]),
    xlabel="Source",
    title="True Amount of Sources in $sample_name Samples",
    clims=(0,1.01)
)

#= True Factor distributions =#
dist_filename = "./data/sundell2022/3Sources from Sundell et al 2022_DISTRIBUTIONS.xlsx"
true_distributions = read_distribution_data(dist_filename);
D_true = to_tensor(true_distributions);

normalize_fiber_sum!(D_true, (2,3));

D_true_processed = mapslices(d -> pre_process(d, n), D_true, dims=1);
F_true = permutedims(D_true_processed,(2,3,1));
#F_true = @view F_true[[1,3,2],:,:]
for (i, F_slice) ∈ enumerate(eachslice(F_true,dims=1))
    p = heatmap(F_slice,
        yticks=(eachindex(F_slice[:,1]), element_list),
        xticks=([1, length(F_slice[1,:])],["0", "1"]),
        xlabel="Normalized Range of Values",
        title = "True Distribution of Elements in Source $i"
        )
    display(p)
end

# Note Ti_temp is off from the factorization...
# but this is because their scales are different.
# Plotting them on the same axis shows they do match:
true_full_scale = true_distributions[3].scale
data_full_scale = distributions[3].scale

true_scale = true_full_scale[1:n÷2:end]
data_scale = data_full_scale[1:n÷2:end]

for i ∈ 1:number_of_factors
    p = plot(true_scale, F_true[i,3,:],
        title = "Learned vs True Ti_temp in source $i",
        xlabel = "value",
        ylabel = "density",
        label = "estimate")
    plot!(data_scale, F[i,3,:], label = "true")
    display(p)
end
