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
    measurements = measurments(D)
    domain_length = length(domains(D)[begin])
    for (name, source) in zip(names(D, 1), eachsource(D))
        full_title = title * "$(dimnames(D)[1]) $name"
        h = Plots.heatmap(
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
    measurment_heatmaps(D::DensityTensor; kw...)

Returns heatmaps for each measurment (lateral slices) of D.

The source name/index for each source will be appended to sourcename.
"""
function measurment_heatmaps(D::DensityTensor; sourcename="", kwargs...) # may clash when looking at a subarray of a DensityTensor
    plots = []
    # No need to normalize since every distribution on the same plot has the same scale
    measurements = measurments(D)
    sources = sources(D)
    for (name, measurement, domain) in zip(measurements(D), eachmeasurment(D), domains(D))
        h = Plots.heatmap(
            domain, sources, array(measurement);
            yflip=true,
            title=title * name,
            kw...
            )
        push!(plots, h)
    end
    return plots
end
