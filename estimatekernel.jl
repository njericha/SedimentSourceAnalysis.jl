#=
Estimate the probability density function for each measurement
=#

using Plots: plot
using KernelDensity
using StatsBase
using OrderedCollections

include("./dataimport.jl")

remove_missing(x) = filter(!ismissing,x) .|> y-> convert(Float64, y)
function inner_percentile!(v, P)
    p_low = (100 - P) / 2
    p_high = 100 - p_low
    a = percentile(v, p_low)
    b = percentile(v, p_high)
    inrange(x) = (a ≤ x ≤ b)
    filter!(inrange, v)
end

function my_bandwidth(v)
    #q25, q75 = quantile(v, [0.25, 0.75])
    #quantile_width = (q75 - q25)/1.34
    #return quantile_width
    return default_bandwidth(v, 1.5)
end

"""
    default_bandwidth(data::AbstractVector{<:Real}, alpha::Float64 = 0.9)

Coppied from KernelDensity since this function is not exported, but I want access
to it so the same bandwidth can be used for different distributions for the
same measurements
"""
function default_bandwidth(data::AbstractVector{<:Real}, alpha::Float64 = 0.9)
    # Determine length of data
    ndata = length(data)
    ndata <= 1 && return alpha

    # Calculate width using variance and IQR
    var_width = std(data)
    q25, q75 = quantile(data, [0.25, 0.75])
    quantile_width = (q75 - q25) / 1.34

    # Deal with edge cases with 0 IQR or variance
    width = min(var_width, quantile_width)
    if width == 0.0
        if var_width == 0.0
            width = 1.0
        else
            width = var_width
        end
    end

    # Set bandwidth using Silverman's rule of thumb
    return alpha * width * ndata^(-0.2)
end

"""
    T, scales, measurements = make_densities(filename::String, n_steps::Integer; P=100)

Creates a tensor of size (n_sinks, n_measurements, n_steps) using data
at filename. Each 3-fiber (there are n_sinks × n_measurements of them) is
a probability density, evenly sampled at n_steps points. The points they are evaluated
at are stored in scales. Note every fiber in the same vertical slice uses the same scale.
The vector measurements stores the names of the measurement for each vertical slice.

P is an optional value between 0 and 100 that filters out each distribution to use only the inner
P percentile range. This can help remove outliers and focus in on where the bulk of the data is.

scales is an optional argument that will use the provided scales instead of the full range of values
for the interpolation points "x".

bandwidths is similar to scales but hold the bandwidth used for each measurement kernal estimation
"""
function make_densities(filename::String, n_steps::Integer; P=100, scales=nothing, bandwidths=nothing) # TODO clean up data types to avoid complicated dictionaries
    sink_data = read_raw_data(filename)
    measurements = collect(keys(sink_data))
    n_measurements = length(measurements)
    n_samples, n_sinks = size(sink_data[measurements[1]]) # TODO more simply extract n_sinks

    # Estimate Kernals Individualy first
    raw_distributions = OrderedDict{String, Vector{UnivariateKDE}}()
    bandwidths_used = zeros(n_measurements)
    for (i, (m, ds)) ∈ enumerate(sink_data)
        kernel_vector = Vector{UnivariateKDE}(undef,n_sinks)
        bandwidth = nothing
        for (j, single_sink_data) ∈ enumerate(eachcol(ds)) # s=1,2,3, but abstracted to allow for more sources
            single_sink_data = remove_missing(single_sink_data)
            inner_percentile!(single_sink_data, P)
            bandwidth = (j==1) ? my_bandwidth(single_sink_data) : bandwidth #ensures the same bandwidth is used for each distribution within the same measurement
            bandwidth = isnothing(bandwidths) ? bandwidth : bandwidths[i]
            kernel_vector[j] = kde(single_sink_data, bandwidth=bandwidth) #kde_lscv(single_sink_data)
        end
        raw_distributions[m] = kernel_vector
        bandwidths_used[i] = bandwidth
    end

    # Standardize Estimations to have the same range within a measurement
    standard_distributions = OrderedDict{String, NamedTuple{(:x,:densities),Tuple{AbstractVector{Float64}, Matrix{Float64}}}}()
    for (i, (m, ds)) ∈ enumerate(raw_distributions)
        # Find extreme values within a measurement
        a = minimum(d -> d.x[begin], ds)
        b = maximum(d -> d.x[end]  , ds)

        # resample each distribution along the larger intervel
        #n_steps = length((ds[1]).x) # use the n_steps given in the input
        x = isnothing(scales) ? range(a, b, n_steps) : scales[i] # use the provided interpolation points if given
        densities = zeros((n_sinks, n_steps))
        for (j, kernel) ∈ enumerate(ds)
            densities[j,:] = pdf(kernel, x) # resample the kernal at the new x values
        end
        standard_distributions[m] = (x=x, densities=densities)
    end

    # Collect data into a single tensor of size (n_sinks, n_measurements, n_steps)
    T = stack((d[2]).densities for d ∈ standard_distributions)::Array{Float64,3}
    T = permutedims(T, [1,3,2])

    # Collect scales for each measurement (the "x" values)
    scales = [(d[2]).x for d ∈ standard_distributions]


    return T, scales, measurements, bandwidths_used
end

"""
    density_draft()

Rough draft of a function to import the raw data and make the density estimates
    as well as plot the estimates.
"""
function density_draft()
    filename_3sources = "./data/sundell2022/3Sources from Sundell et al 2022.xlsx"
    filename_20sinks = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"

    # Start by observing the source data
    true_data = read_raw_data(filename_3sources)
    measurements = collect(keys(true_data))

    # Plot distribution for each measurement and source
    # These are the ground truth distributions
    for m ∈ measurements
        for (i, single_source_data) ∈ enumerate(eachcol(true_data[m])) # s=1,2,3, but abstracted to allow for more sources
            single_source_data = remove_missing(single_source_data)
            kernel = kde(single_source_data)
            p = plot(kernel.x, kernel.density,
                title="$m Probability Density Estimate for Source $i",
                label="pdf estimate",
                xlabel="value",
                ylabel="density")
            scatter!(single_source_data, 0 .* single_source_data,
                label="data points")
            display(p)
        end
    end

    # Generate distribution for each measurement in the 20 sinks/rocks
    sink_data = read_raw_data(filename_20sinks)
    measurements = collect(keys(sink_data))

    n_samples, n_sinks = size(sink_data[measurements[1]])
    n_measurements = length(measurements)
    raw_distributions = Dict{String, Vector{UnivariateKDE}}()

    showplots = false

    for m ∈ measurements
        kernel_vector = Vector{UnivariateKDE}(undef,n_sinks)
        for (i, single_sink_data) ∈ enumerate(eachcol(sink_data[m])) # s=1,2,3, but abstracted to allow for more sources
            single_sink_data = remove_missing(single_sink_data)
            kernel = kde(single_sink_data)
            p = plot(kernel.x, kernel.density,
                title="$m Probability Density Estimate for Sink $i",
                label="pdf estimate",
                xlabel="value",
                ylabel="density")
            scatter!(single_sink_data, 0 .* single_sink_data,
                label="data points")
            showplots ? display(p) : nothing
            kernel_vector[i] = kernel
        end
        raw_distributions[m] = kernel_vector
    end

    # Standardize distribution ranges for the same measurement between sinks
    standard_distributions = Dict{String, NamedTuple{(:x,:densities),Tuple{AbstractVector{Float64}, Matrix{Float64}}}}()

    for (measurement, distributions) ∈ collect(raw_distributions)
        # find minimum and maximum "x" value across the sinks
        a = minimum(d -> d.x[begin], distributions)
        b = maximum(d -> d.x[end]  , distributions)

        # resample each distribution along the larger intervel
        n_steps = length((distributions[1]).x)
        x = range(a, b, n_steps)
        densities = zeros((n_sinks, n_steps))
        for (i, kernel) ∈ enumerate(distributions)
            density = pdf(kernel, x) # resample the kernal at the new x values
            densities[i,:] = density

            p = plot(x, density,
                title="Standardized $measurement Density for Sink $i",
                label="pdf estimate",
                xlabel="value",
                ylabel="density",
            )
            showplots ? display(p) : nothing
        end
        standard_distributions[measurement] = (x=x, densities=densities)
    end

    # Collect all distributions into a single 3-order tensor

    # function to_tensor(D::Dict{String, NamedTuple{(:x,:densities),Tuple{AbstractVector{Float64}, Matrix{Float64}}}})
    #     T = zeros(())
    # end
    # densities_stack = Matrix[]
    # for d ∈ standard_distributions
    #     push!(densities_stack, (d[2]).densities)
    # end
    # n_steps = length((standard_distributions[measurements[1]]).x)
    # T = zeros((n_sinks, n_measurements, n_steps))
    T = stack((d[2]).densities for d ∈ standard_distributions)
    T = permutedims(T, [1,3,2])
    scales = [(d[2]).x for d ∈ standard_distributions]
    return T, scales
end
