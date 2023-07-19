function inner_percentile!(v, P)
    p_low = (100 - P) / 2
    p_high = 100 - p_low
    a = percentile(v, p_low)
    b = percentile(v, p_high)
    inrange(x) = (a ≤ x ≤ b)
    filter!(inrange, v)
end

global DEFAULT_ALPHA = 0.9::Real

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
    make_distributions(s::Sink; inner_percentile::Integer=100, bandwidths::AbstractVector{<:Real})

Estimates the distributions for each measurment in a Sink
returns the UnivariateKDE as well as samples of the kernel

# Parameters
- `inner_percentile`: value between 0 and 100 that filters out each distribution
and uses the inner percentile range. This can help remove outliers and focus in on where
the bulk of the data is.
- `bandwidths`: list of bandwidths used for each measurement's density estimation

# Returns
- `density_estimates`
"""
function make_distributions(
    s::Sink;
    inner_percentile::Integer=100,
    bandwidths::AbstractVector{<:Real},
    )
    # Argument Handeling
    # Check input is in the correct range
    (0 < inner_percentile <= 100) ||
        ArgumentError("inner_percentile must be between 0 and 100, got $inner_percentile")

    # Loop setup
    eachmeasurment = eachmeasurment(s)
    n_measurments = length(eachmeasurment)
    density_estimates = Vector{UnivariateKDE}(undef, n_measurments)

    for (i, (measurment_values, b)) in enumerate(zip(eachmeasurment, bandwidths))
        # Estimate density based on the inner precentile to ignore outliers
        inner_percentile!(measurment_values, inner_percentile)
        KDEs[i] = kde(measurment_values, bandwidth=b)
    end

    return density_estimates
end

function make_distributions(s::Sink; inner_percentile::Integer=100, alpha=DEFAULT_ALPHA)
    bandwidths = default_bandwidth.(eachmeasurment(s), alpha)
    return (make_distributions(
        s::Sink;
        inner_percentile=inner_percentile,
        bandwidths=bandwidths,
        ), bandwidths)
end

const DEFAULT_N_SAMPLES = 64::Integer

"""
    standardize_KDEs(KDEs::AbstractVector{UnivariateKDE}; n_samples=DEFAULT_N_SAMPLES,)

Resample the distributions within each sink so that like-measurments use the same scale
"""
function standardize_KDEs(KDEs::AbstractVector{UnivariateKDE}; n_samples=DEFAULT_N_SAMPLES,)
    a = minimum(d -> d.x[begin], KDEs) # smallest left endpoint
    b = maximum(d -> d.x[end]  , KDEs) # biggest right endpoint

    x_new = range(a, b, length=n_samples) # make the (larger) x-values range
    KDEs_new = pdf.(KDEs, x_new) # resample the densities on the new range

    return KDEs_new, x_new
end
