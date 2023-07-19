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

global DEFAULT_ALPHA = 0.9

function set_alpha(alpha)
    DEFAULT_ALPHA = alpha
end

"""
    default_bandwidth(data::AbstractVector{<:Real}, alpha::Float64 = 0.9)

Coppied from KernelDensity since this function is not exported. I want access
to it so that the same bandwidth can be used for different distributions for the
same measurments.
"""
function default_bandwidth(data::AbstractVector{<: Real}, alpha::Float64 = DEFAULT_ALPHA)
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
    make_densities(filename::String, n_steps::Integer; P=100)

Creates a tensor of size (n_sinks, n_measurements, n_steps) using data
at filename. Each 3-fibre (there are n_sinks × n_measurements of them) is
a probability density, evenly sampled at n_steps points. The points they are evaluated
at are stored in scales. Note every fibre in the same vertical slice uses the same scale.
The vector measurements stores the names of the measurment for each vertical slice.

P is an optional value between 0 and 100 that filters out each distribution to use only the inner
P percentile range. This can help remove outliers and focus in on where the bulk of the data is.

scales is an optional argument that will use the provided scales instead of the full range of values
for the interpolation points "x".

bandwidths is similar to scales but hold the bandwidth used for each measurement kernel estimation

# Returns
T, scales, measurements, bandwidths_used
"""
function make_distributions(d::Sink, n_steps::Integer; P=100, scales=nothing, bandwidths=nothing) # TODO clean up data types to avoid complicated dictionaries
    sink_data = d
    measurements = collect(keys(sink_data))
    n_measurements = length(measurements)
    n_samples, n_sinks = size(sink_data[measurements[1]]) # TODO more simply extract n_sinks

    # Estimate Kernels Individualy first
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
            densities[j,:] = pdf(kernel, x) # resample the kernel at the new x values
        end
        standard_distributions[m] = (x=x, densities=densities)
    end

    # Collect data into a single tensor of size (n_sinks, n_measurements, n_steps)
    T = stack((d[2]).densities for d ∈ standard_distributions)::Array{Float64,3}
    T = permutedims(T, [1,3,2])

    # Collect scales for each measurment (the "x" values)
    scales = [(d[2]).x for d ∈ standard_distributions]


    return T, scales, measurements, bandwidths_used
end

"""
    make_distributions(s::Sink)

estimates the distributions for each measurment in a DataSink
returns the UnivariateKDE as well as samples of the kernel
"""
function make_distributions(s::DataSink; n_samples=64, inner_percentile=100, bandwidths=nothing)
    # Check input is in the correct range
    (0 < inner_percentile <= 100) ||
        ArgumentError("inner_percentile must be between 0 and 100, got $inner_percentile")

    # TODO

end

function make_distribution(v, n_samples, inner_percentile)
end

"""
    standardize!(sinks::AbstractVector{Sink})
    standardize!(sink1, sink2, ...)

Resample the distributions within each sink so that like-measurments use the same scale
"""
function standardize!(sinks::AbstractVector{DistributionSink})
    # Get the scales for each sink
    # For each measurment
        # Find the outer ranges of its scale (the x values)
        # resample the values on this (larger) range
        # use the values to sample the distributions
        # Save the newly sampled distribution to that sink
end

standardize!(sinks::AbstractVector{DistributionSink}...) = standardize!(sinks)
