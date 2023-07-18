#=
    sortgrain.jl

Functions for taking a single grain data and finding which
source it most likely belongs to.
=#
using OrderedCollections

# Transform data to a more useful format by grain
# Or have a function look up the vector of measurments for the grain
"""
    get_grain(raw_data, (grain_i, sink_j))

obtains grain_i from sink_j in the raw_data and returns a vector
of the measurments of that grain
"""
function get_grain(raw_data, (grain_i, sink_j))
    n_measurements = length(raw_data)
    grain_vec = zeros(n_measurements)
    for (k, (_, measurement)) ∈ enumerate(raw_data)
        grain_vec[k] = measurement[grain_i, sink_j]
    end
    return grain_vec
end

# Compare with the distribution for a factor
"""
    find_subinterval(scale, grain_value)

Find which subinterval the data lies in for each measurement

```jdocstest
julia> test_scale = 1:10
1:10

julia> find_subinterval(test_scale, 4.8)
5
```
"""
function find_subinterval(scale, grain_value)
    val, k = findmin(x->abs(x-grain_value), scale)
    return k
end



"""
    estimate_prob(factor_slice, scales, grain_vec)

Estimates the probability of observing grain_vec from the single factor
This assumes the distributions are normalized in the sence that
they carry equal weight e.i. the row_sums are 1.

Note the distributions before rescaling have a row_sum equal to 1/stepsize
with stepsize being the intervel width used in scale.
"""
function estimate_prob(factor_slice, scales, grain_vec)
    measurment_probabilities = zeros(n_measurements)
    for (j, (grain_value, density, scale)) ∈ enumerate(zip(grain_vec, eachrow(factor_slice), scales))
        subinterval_k = find_subinterval(scale, grain_value)
        density_k = density[subinterval_k] # Get the probability for observing a measurment on that intervel
        measurment_probabilities[j] = density_k
    end
    ## Multiply the probabilities to get the overall probability of
    ## observing the data in that volume
    return prod(measurment_probabilities)
end

"""
    (prob, factor_i) = estimate_which_factor(F, scales, grain_vec)

Returns the likelyhood and factor number of the mostly likely factor
the grain vector came from
"""
function estimate_which_factor(F, scales, grain_vec)
    n_factors, _, _ = size(F)
    likelyhoods = zeros(n_factors)
    for (i, factor_slice) ∈ enumerate(eachslice(F, dims=1)) #factor_slice = F[1, :, :]
        prob = estimate_prob(factor_slice, scales, grain_vec)
        likelyhoods[i] = prob
    end
    return findmax(likelyhoods), likelyhoods #find the one with the highest probability
end

"""
    estimate_factors_and_probs(raw_data, sink_j)

Takes the raw data as an ordered dictionary and an integer corresponding
to the jth sink and returns the guess for which factor each grain in that sink
belongs to, and the probability estimate.
"""
function estimate_factors_and_probs(F, scales, raw_data, sink_j)
    n_grains = 75 #, n_sinks = size(iterate(raw_data)[1][2])[1]
    factor_estimates = zeros(Int, n_grains)
    prob_estimates = zeros(n_grains)
    for i ∈ 1:n_grains # TODO use getindex
        grain_vec = get_grain(raw_data, (i, sink_j))
        (prob, factor_i), _ = estimate_which_factor(F, scales, grain_vec)
        factor_estimates[i] = factor_i
        prob_estimates[i] = prob
    end
    return factor_estimates, prob_estimates
end

# Shorthand Helpers
estimate_factor_probs(grain_i,sink_j) = (estimate_which_factor(F, scales, get_grain(raw_data, (grain_i,sink_j))))[2]
factor_probs(grain_i,sink_j) = (estimate_which_factor(F_true, scales, get_grain(raw_data, (grain_i,sink_j))))[2]

#= Start of Script Portion of File =#

# Load the learned coefficients C and factor distributions F
using JLD2
using UnPack

jldopen("estimated_factors.jld2", "r") do file
    @unpack C, F, scales, names = file
end

# Import the raw grain data
include("./dataimport.jl")
filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
raw_data = read_raw_data(filename)

n_sinks = 20
n_sources = 3

# Learned distributions
for sink_j ∈ 1:n_sinks
    # Estimate the factors and probs for the first sink
    factor_estimates, prob_estimates = estimate_factors_and_probs(
        F, scales, raw_data, sink_j
    )

    # Plot the factor estimates, where brighter colors indicate a higher probability
    p = scatter(factor_estimates;
        zcolor=log10.(prob_estimates),
        clims=(-8,-6),
        yticks = 1:n_sources,
        ylabel = "source",
        xlabel = "grain index",
        title = "Source Estimate and probabilities for Sink $sink_j",
        zlabel = "10^_"
    )
    display(p)
end

# True distributions
dist_filename = "./data/sundell2022/3Sources from Sundell et al 2022.xlsx"
F_true, scales_true, names_true, bandwidths_true = make_densities(dist_filename, 2^5, P=100, scales=scales, bandwidths=bandwidths);

Δxs = [scale[2]-scale[1] for scale ∈ scales] # TODO replace this with scale_slices!() once the function is moved and avalible
ds_true = eachslice(F_true, dims=2)
ds_true .*= Δxs#slice_norm

for sink_j ∈ 1:n_sinks
    # Estimate the factors and probs for the first sink
    factor_estimates, prob_estimates = estimate_factors_and_probs(
        F_true, scales, raw_data, sink_j
    )

    # Plot the factor estimates, where brighter colors indicate a higher probability
    p = scatter(factor_estimates;
        zcolor=log10.(prob_estimates),
        clims=(-8,-6),
        yticks = 1:n_sources,
        ylabel = "source",
        xlabel = "grain index",
        title = "True Source Estimate and Probabilities for Sink $sink_j",
        zlabel = "10^_"
    )
    display(p)
end

# Ground truth grain sources
source_filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
source_proportions = XLSX.readdata(source_filename,"Source proportions","B2:D21")
source_proportions = convert.(Int, source_proportions)
"""
    make_source_vector(integers::AbstractVector{Integer})

Given a list of integers, make a new list containing
integers[1] many 1s, integers[2] many 2s, etc.

Example
-------
make_source_vector([3,4,2]) == [1 1 1 2 2 2 2 3 3]
"""
function make_source_vector(integers::AbstractVector{T}) where T  <: Integer
    v = Int[]
    for (i, n) ∈ enumerate(integers)
        v_i = fill(i, n)
        append!(v, v_i)
    end
    return v
end

for (i, proportion) ∈ enumerate(eachrow(source_proportions))
    true_factors = make_source_vector(proportion)

    # Plot the factor estimates, where brighter colors indicate a higher probability
    p = scatter(true_factors;
    yticks = 1:n_sources,
    ylabel = "source",
    xlabel = "grain index",
    title = "True Sources for Sink $i"
    )
    display(p)
end
