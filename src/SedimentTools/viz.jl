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
    measurements = measurements(D)
    domain_length = length(domains(D)[begin])
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
function measurement_heatmaps(D::DensityTensor; kwargs...) # may clash when looking at a subarray of a DensityTensor
    plots = []
    # No need to normalize since every distribution on the same plot has the same scale
    measurements = measurements(D)
    sources = sources(D)
    for (name, measurement, domain) in zip(measurements(D), eachmeasurement(D), domains(D))
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
    measurement_heatmaps(D::DensityTensor, measurement::String; kw...)

Returns one plot will all distributions for a given measurement (lateral slice) of D.
"""
function measurement_distributions(D::DensityTensor, measurement::String; kwargs...) # may clash when looking at a subarray of a DensityTensor
    # No need to normalize since every distribution on the same plot has the same scale
    domain = domain(D, measurement)
    p = plot()
    for (source_name, density) in zip(sources(D), eachdensity(D, measurement))
        plot!(domain, density; label=source_name, kw...)
    end
    return p
end
