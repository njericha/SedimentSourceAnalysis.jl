"""
    Plots.heatmap(M::NamedMatrix; kwargs...)

Automatically grabs the names in axis to label the x and y plot axis
"""
function Plots.heatmap(M::NamedMatrix; kwargs...) # may clash when looking at a subarray of a DensityTensor
    return Plots.heatmap(names(M,2), names(M,1), M.array; yflip=true, kwargs...)
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
            collect(source);
            yticks=(eachindex(measurements), measurements),
            xticks=([1, domain_length],["min", "max"]),
            xlabel="Typical Range of Values",
            yflip=true,
            title=full_title,
            kwargs...
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
    sources = getsourcenames(D)
    sourcename = getsourcename(D)
    for (name, measurement, domain) in zip(getmeasurements(D), eachmeasurement(D), getdomains(D))
        h = heatmap(
            domain, sources, array(measurement);
            yticks=(eachindex(sources), sources),
            yflip=true,
            title=title * name,
            ylabel=sourcename,
            kwargs...
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

"""
    plot_source_index(
    indexes::AbstractVector{<:Integer},
    loglikelihood_ratios::AbstractVector{<:Real};
    kwargs...
    )

Returns one scatter plot of dots (eachindex(indexes), indexes) with brighter colours
corresponding to higher loglikelihood_ratios.
"""
function plot_source_index(
    indexes::AbstractVector{<:Integer},
    loglikelihood_ratios::AbstractVector{<:Real};
    kwargs...
    )
    p = scatter(
            indexes;
            marker_z=loglikelihood_ratios,
            clabel="log likelihood ratio",
            colorbar_ticks=(0:0.5:1, ["0", "0.5", ">1.0"]),
            clim=(0,1),
            xlabel="grain index",
            ylabel="source",
            yticks=sort([Set(indexes)...]), # sort the unique indexes that appear
            label="",
            kwargs...
    )
    return p
end

"""
    plot_convergence(rel_errors, norm_grad, dist_Ncone)

Returns 3 separate plots for the three convergence metrics on a log10-y scale.
"""
function plot_convergence(rel_errors, norm_grad, dist_Ncone)
    return plot_convergence((rel_errors, norm_grad, dist_Ncone))
end

"""
    plot_convergence((rel_errors, norm_grad, dist_Ncone))
"""
function plot_convergence(all_series::NTuple{3, AbstractVector{<:Real}})
    titles = (
        "Relative Error Between Y and C*F Convergence",
        "Norm of Full Gradient Convergence",
        "Distance Between Negative Gradient and Normal Cone",
    )

    ylabels = (
        "Relative Error",
        "Norm of Gradient",
        "Distance to Normal Cone",
    )

    plots = []
    for (series, title, ylabel) in zip(all_series, titles, ylabels)
        p = plot(
            series;
            ylabel,
            title,
            xlabel="iteration #",
            yscale=:log10,
            legend=false,
        )
        push!(plots, p)
    end

    return plots
end
