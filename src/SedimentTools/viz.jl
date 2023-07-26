"""
    Plots.heatmap(M::NamedMatrix; kwargs...)

Automatically grabs the names in axis to label the x and y plot axis
"""
function Plots.heatmap(M::NamedMatrix; kwargs...) # may clash when looking at a subarray of a DensityTensor
    return Plots.heatmap(names(M,2), names(M,1), M.array; yflip=true, kw...)
end

"""
    source_heatmaps(D::DensityTensor; title="", kw...)

Returns heatmaps for each source (horizontal slices) of D.

The source name/index for each source will be appended to title.
"""
function source_heatmaps(D::DensityTensor; title="", kwargs...) # may clash when looking at a subarray of a DensityTensor
    plots = []
    D = normalize_density_sums(D)
    measurements = getmeasurements(D)
    domain_length = length(getdomains(D)[begin])
    for (name, source) in zip(names(D, 1), eachsource(D))
        full_title = title * "$(dimnames(D)[1]) $name"
        h = heatmap(
            array(source);
            yticks=(eachindex(measurements), measurements),
            xticks=([1, domain_length],["min", "max"]),
            xlabel="Typical Range of Values",
            yflip=true,
            title=title * name,
            kw...
            )
        push!(plots, h)
    end
    return plots
end

"""
    measurement_heatmaps(D::DensityTensor; kw...)

Returns heatmaps for each measurement (lateral slices) of D.
"""
function measurement_heatmaps(D::DensityTensor; title="", kwargs...) # may clash when looking at a subarray of a DensityTensor
    plots = []
    # No need to normalize since every distribution on the same plot has the same scale
    measurements = getmeasurements(D)
    sources = getsourcenames(D)
    for (name, measurement, domain) in zip(getmeasurements(D), eachmeasurement(D), getdomains(D))
        h = heatmap(
            domain, sources, array(measurement);
            yflip=true,
            title=title * name,
            kw...
            )
        push!(plots, h)
    end
    return plots
end

"""
    distributions_plot(D::DensityTensor, measurement::String; kw...)

Returns one plot will all distributions for a given measurement (lateral slice) of D.
"""
function plot_densities(D::DensityTensor, measurement::String; kwargs...) # may clash when looking at a subarray of a DensityTensor
    # No need to normalize since every distribution on the same plot has the same scale
    domain = getdomain(D, measurement)
    p = plot(ylabel="density")
    for (source_name, density) in zip(getsourcenames(D), eachdensity(D, measurement))
        plot!(domain, density; label=source_name, kwargs...)
    end
    return p
end
