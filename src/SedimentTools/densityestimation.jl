"""
    make_densities(s::Sink; kwargs...)
    make_densities(s::Sink, domains::AbstractVector{<:AbstractVector}; kwargs...)

Estimates the densities for each measurement in a Sink.

When given domains, a list where each entry is a domain for a different measurement,
resample the kernel on this domain.

# Parameters
- `bandwidths::AbstractVector{<:Real}`: list of bandwidths used for each measurement's
density estimation
- `inner_percentile::Integer=100`: value between 0 and 100 that filters out each measurement
by using the inner percentile range. This can help remove outliers and focus in on where the
bulk of the data is.

# Returns
- `density_estimates::Vector{UnivariateKDE}`
"""
function MatrixTensorFactor.make_densities(
    s::Sink;
    inner_percentile::Integer=100,
    bandwidths::AbstractVector{<:Real}=default_bandwidth.(
        collect(eachmeasurement(s)),DEFAULT_ALPHA,inner_percentile),
    )
    # Argument Handeling: check inner_percentile is a percentile
    (0 < inner_percentile <= 100) ||
        ArgumentError("inner_percentile must be between 0 and 100, got $inner_percentile")

    # Loop setup
    _eachmeasurement = eachmeasurement(s)
    n_measurements = length(_eachmeasurement)
    density_estimates = Vector{UnivariateKDE}(undef, n_measurements)

    for (i, (measurement_values, b)) in enumerate(zip(_eachmeasurement, bandwidths))
        # Estimate density based on the inner precentile to ignore outliers
        measurement_values = filter_inner_percentile(measurement_values, inner_percentile)
        density_estimates[i] = kde(measurement_values, bandwidth=b)
    end

    return density_estimates
end

function MatrixTensorFactor.make_densities(
    sinks::Sink,
    domains::AbstractVector{<:AbstractVector};
    kwargs...
    )
    KDEs = MatrixTensorFactor.make_densities(sinks; kwargs...)
    KDEs_new = pdf.(KDEs, domains)
    return KDEs_new
end
