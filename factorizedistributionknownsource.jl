#=
This script uses the distributions made by KernelDensity.jl
=#

using Plots
using Statistics
using StatsPlots

using Printf # for formating floats as strings
using JLD2, UnPack # for saving and loading data to a file


include("./coorddecent.jl")
include("./estimatekernel.jl")

# Import data
#include("./dataimport.jl")
filename = "./sundelldata/20sinks from 3Sources from Sundell et al 2022.xlsx"
T, scales, names, bandwidths = make_densities(filename, 2^5, P=95);

# Change measurment name "ages" to "Age" for consistantcy with the measurements in the sources file
age_idx = findall(s -> s == "ages", names)
names[age_idx] .= "Age"


# Look at distributions
for (j, (scale, name)) ∈ enumerate(zip(scales,names))
    p = plot(scale,eachrow(T[:,j,:]),
        title="$name Distribution for Each Sink",
        yaxis="density",
        xaxis= "value")
    display(p)
end

# Ensure every distribution has area 1 on average by rescaling by the stepsize used to
# evalutate the distribution estimate
function scale_slices!(T, scales, dims)
    Δx = [scale[2]-scale[1] for scale ∈ scales]
    d = eachslice(T, dims=dims)
    d .*= Δx
end

scale_slices!(T, scales, 2)
isapprox_one(x) = isapprox(x, 1.0; atol=1e-4)
@assert all(isapprox_one.(sum.(eachslice(T,dims=(1,2)))))

# Perform Decomposition
normalize_each_update = true
number_of_factors = 3
C, F, error, norm_grad, dist_Ncone = coorddecent(T, number_of_factors,
    maxiter=10000, tol=1e-4, normalize_each_update=normalize_each_update,
    plot_progress = 0); #1e-5 for norm_grad convergence #1e-4 for Normal cone convergence with normalizing each update

# Rescale learned sources, if normalization was not done at every step
if !normalize_each_update
    fiber_sums = sum.(eachslice(F,dims=(1,2)))
    avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
    F = avg_factor_sums^(-1) * F # TODO make more accurate scaling
    C = C * avg_factor_sums
end

# Import ground truth matricies C and F
source_filename = "./sundelldata/20sinks from 3Sources from Sundell et al 2022.xlsx"
source_amounts = XLSX.readdata(source_filename,"Source proportions","B2:D21")
C_true = source_amounts / 75 # 75 Grains in each sink

dist_filename = "./sundelldata/3Sources from Sundell et al 2022.xlsx"
F_true, scales_true, names_true, bandwidths_true = make_densities(dist_filename, 2^5, P=100, scales=scales, bandwidths=bandwidths);
scale_slices!(F_true, scales, 2)

# Check that all distributions are in the same order with the same scaling and bandwidths
@assert scales_true == scales
@assert names_true == names
@assert bandwidths_true == bandwidths

"""
    permute_factors!(C, F, C_true, F_true)

Permute factors in C and F to match the ground truth C_true and F_true
"""
function permute_factors!(C, F, C_true, F_true; double_check=false)
    # make list to store which true source the columns of C match
    n_factors = size(C_true)[2]
    true_ordering = zeros(Integer, n_factors)

    # loop over every column of C and find the best matching column of C_true
    for (i, c_true) ∈ enumerate(eachcol(C_true))
        _, i_true = findmin(c -> norm(c - c_true), eachcol(C))
        true_ordering[i] = i_true
    end
    @assert allunique(true_ordering)

    if double_check #repeat for F
        true_ordering2 = zeros(Integer, n_factors)
        for (i, slice_true) ∈ enumerate(eachslice(F_true, dims=1))
            _, i_true = findmin(s -> norm(s - slice_true), eachslice(F, dims=1))
            true_ordering2[i] = i_true
        end
        @assert true_ordering == true_ordering2
    end

    # Swap columns of C and horizontal slices of F to the new ordering
    C .= @view C[:,true_ordering]
    F .= @view F[true_ordering,:,:]
    return true_ordering
end
@time true_ordering = permute_factors!(C, F, C_true, F_true)

# Save C and F to a file
jldopen("estimated_factors.jld2", "w") do file
    @pack! file = C, F, scales, names, bandwidths
end

# can open the file like so:
# jldopen("estimated_factors.jld2", "r") do file
#    @unpack C, F, scales, names, bandwidths = file
# end

#### Plotting ####

# Plot gradient convergence
plot(norm_grad ./ norm_grad[1],
    yaxis=:log10,
    ylabel="2-norm of Full Gradient",
    title="Convergence of Gradient",
    xlabel="iteration #")

plot(dist_Ncone,
    yaxis=:log10,
    ylabel="distance",
    title="Distance of -gradient to Normal Cone of Constraint",
    xlabel="iteration #")

# Plots the resulting coefficients
sample_name = "sink"
heatmap(C,
    yticks=(eachindex(C[:,1]), ["$sample_name $i" for i ∈ eachindex(C[:,1])]),
    xticks=eachindex(C[1,:]),
    xlabel="Factor",
    title="Learned Amount of Factors in $sample_name Samples",
    clims=(0,1.01)
)

# For each element, compare distributions between factors
for (j, (scale, name)) ∈ enumerate(zip(scales,names))
    F_slice = F[:,j,:]
    #name = short_name(name)
    #xlabel = (lowercase(name) == "age") ?  "Ma" : "ppm"
    xlabel = "value"
    p = heatmap(scale, (eachindex(F_slice[:,1])), F_slice,
        xlabel = xlabel,
        ylabel = "Factor",
        title = "Distribution of $name in Factors",
        yticks = eachindex(F_slice[1,:])
    )
    display(p)
end

# Show distribution of different elements within a factor
for (i, F_slice) ∈ enumerate(eachslice(F,dims=1))
    p = heatmap(F_slice,
        yticks=(eachindex(F_slice[:,1]), names),
        xticks=([1, length(F_slice[1,:])],["0", "1"]),
        xlabel="Normalized Range of Values",
        title = "Learned Distribution of Elements in Factor $i"
        )
    display(p)
end

#= True Source proportions =#

heatmap(C_true,
    yticks=(eachindex(C_true[:,1]), ["$sample_name $i" for i ∈ eachindex(C_true[:,1])]),
    xticks=eachindex(C_true[1,:]),
    xlabel="Source",
    title="True Amount of Sources in $sample_name Samples",
    clims=(0,1.01)
)

#= True Factor distributions =#

for (i, F_slice) ∈ enumerate(eachslice(F_true,dims=1))
    p = heatmap(F_slice,
        yticks=(eachindex(F_slice[:,1]), names),
        xticks=([1, length(F_slice[1,:])],["0", "1"]),
        xlabel="Normalized Range of Values",
        title = "True Distribution of Elements in Source $i"
        )
    display(p)
end