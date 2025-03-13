#=
Use data from /data/sundell2022

Try using the median of bandwidths for each feature to do the KDEs
rather than taking the bandwidth of the first sink
=#

using XLSX
using NamedArrays
using Plots
using Printf
using Random
using Statistics: mean, median
using LinearAlgebra: norm
using Logging; disable_logging(Warn)

using BlockTensorDecomposition
using SedimentSourceAnalysis

# Plot settings
Plots.resetfontsizes(); Plots.scalefontsizes(1.5)
plotfont="Computer Modern"
plotfontsize=13
default(legendfontsize=plotfontsize, plot_titlefontsize=plotfontsize+1, titlefont=plotfontsize, legendtitlefontsize=plotfontsize, fontfamily=plotfont, legendfont=plotfont)

# set random seed for repeatability
# does not randomize data, but the initialization for nnmtf
Random.seed!(314159265)

# Import data from excel file
filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
sinks = read_raw_data(filename)::Vector{Sink}

## Look at a grain
sink1 = sinks[begin]
grain1 = sink1[begin]
println("Grain 1 in Sink 1")
display(grain1)

## Here we see each measurement displayed as well as the sampled value
## To get the names of each measurement, use getmeasurements(.)
@show getmeasurements(grain1)

# Estimate the densities of each sink

## Select the bandwidth for the estimation
## Uses Silverman's rule of thumb
sink1 = sinks[begin]
inner_percentile = 95 # Filter outliers; ignore values outside the inner percentile
alpha_ = 0.9 # bandwidth "alpha" smooths the density estimate, 0.9 is the default
             # this can denoise estimation

bandwidth_matrix = zeros(length(sinks), length(grain1))
for (i, sink) in enumerate(sinks)
    bandwidth_matrix[i, :] = default_bandwidth.(collect(eachmeasurement(sink)), alpha_, inner_percentile)
end

bandwidths = dropdims(median(bandwidth_matrix, dims=1), dims=1) # column-wise median

## Note getmeasurements() gets the *names* of each measurement,
## whereas eachmeasurement() is an iterator for the *values* for each measurement

## Obtain the raw densities estimates
## The same measurement could (and likely!) have different supports for different sinks...
raw_densities = make_densities.(sinks; bandwidths, inner_percentile)

## ...so we standardize them by resampling the densities on the same domain,
## which is the union of all intervels.
densities, domains = standardize_KDEs(raw_densities) #Made it to here so far

## We now have a list of densities and domain xs they are sampled at.
## Here, densities[i] are the densities for the ith sink (sinks[i]),
## and xs[j] is the domain for the jth measurement (densities[i][j]).
## Importantly, this is the same domain regardless of the sink!

# Package the densities into a single order-3 tensor
densitytensor = DensityTensor(densities, domains, sinks);
setsourcename!(densitytensor, "sink");

# Visualize the data in the tensor by plotting the densities for the first measurement
measurement_names = getmeasurements(densitytensor) # ["Age", ...]

function plot_feature(name, x_label)
    p = plot_densities(densitytensor, name);
    plot!(xlabel=x_label,
        ylabel="probability density",
        legendtitle = "Sink",
        legend_columns=2,
        )
    display(p)

    # Vizualize the first sink's KDE for the first measurement (Age)
    raw_data_ages_sink1 = [grain[name] for grain in sink1]
    KDE_ages_sink1 = densitytensor[1, name, :] #sink 1
    p = histogram(raw_data_ages_sink1;
        normalize=true,
        #bins=range(0, 200, length=6),
        label="histogram",
        alpha=0.25,
        color=:black,
        )
    plot!(getdomain(densitytensor, name), KDE_ages_sink1;
        label="KDE",
        color=:blue,
        linewidth=5,
        alpha=1
        )
    scatter!(raw_data_ages_sink1, zeros(length(raw_data_ages_sink1));
        marker=:vline,
        markersize=15,
        label="grain sample",
        color=:black,
        xlabel=x_label,
        ylabel="probability density",
        )
    display(p)

    # Vizualize the first sink's KDE and sample points
    p = plot(getdomain(densitytensor, name), KDE_ages_sink1;
    label="continuous KDE",
    color=:blue,
    xlabel=x_label,
    ylabel="probability density"
    )
    scatter!(getdomain(densitytensor, name), KDE_ages_sink1;
    marker=:circ,
    markercolor=:blue,
    label="KDE discretization")
    display(p)

    # Coarse Vizualize the first sink's KDE and sample points
    p = plot(getdomain(densitytensor, name), KDE_ages_sink1;
    label="continuous KDE",
    color=:blue,
    xlabel=x_label,
    ylabel="probability density"
    )
    scatter!(getdomain(densitytensor, name)[1:4:end], KDE_ages_sink1[1:4:end];
    marker=:circ,
    markercolor=:blue,
    label="KDE discretization")
    display(p)

    # Coarse Vizualize the first sink's KDE and sample points
    p = plot(getdomain(densitytensor, name), KDE_ages_sink1;
    label="continuous KDE",
    color=:blue,
    xlabel=x_label,
    ylabel="probability density"
    )

    display(p)

    n_boxes = 14
    step = length(KDE_ages_sink1) ÷ n_boxes
    y_coarse=KDE_ages_sink1[1:step:end]
    x_coarse=getdomain(densitytensor, name)[1:step:end]

    # Get sample points
    scatter!(x_coarse, y_coarse, label="discretized KDE")

    display(p)

    function make_intervals(xs)
        return [[xs[i],xs[i+1]] for i in 1:length(xs)-1]
    end

    rectangle(x1, y1, x2, y2) = Shape([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])

    # Plot Riemann Sums
    intervals = make_intervals(x_coarse)
    for (i, X) in enumerate(intervals)
        i == 1 ? continue : nothing
        Δx = (X[2]-X[1])/2
        box = rectangle(X[1]-Δx,0,X[2]-Δx,y_coarse[i])
        i == 2 ? label = "probability values" : label = nothing
        plot!(box;c=:blue,label,alpha=.3)
    end

    display(p)
end


plot_feature("Eu_anomaly", "Eu Anomaly")

plot_feature("Ti_temp", "Ti Temperature")
