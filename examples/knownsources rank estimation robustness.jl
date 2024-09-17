#=
Use data from /data/sundell2022

Perform the rank estimation using different
1) bandwidths
2) number of grains in each sink
=#

using XLSX
using NamedArrays
using Plots
using Printf
using Random
using Statistics: mean, median
using LinearAlgebra: norm
using Logging; disable_logging(Warn)

using MatrixTensorFactor
using SedimentAnalysis

# Plot settings
Plots.resetfontsizes(); Plots.scalefontsizes(1.5)
plotfont="Computer Modern"
plotfontsize=13
default(legendfontsize=plotfontsize, plot_titlefontsize=plotfontsize+1, titlefont=plotfontsize, legendtitlefontsize=plotfontsize, fontfamily=plotfont, legendfont=plotfont)

# set random seed for repeatability
# does not randomize data, but the initialization for nnmtf
Random.seed!(314159265)

# Tuple extracting
compute_final_loss_curve((scale_bandwidth, grains_per_sink)) = compute_final_loss_curve(scale_bandwidth, grains_per_sink)

function compute_final_loss_curve(scale_bandwidth, proportion_of_grains, n_features=7, n_sinks=20, discretization_size=64) # 7 features and 20 sinks in the dataset
    skip_grain = (1 / proportion_of_grains) |> round |> Int

        # Import data from excel file
    filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
    sinks = read_raw_data(filename)::Vector{Sink}

    sinks = [sink[begin:skip_grain:end] for sink in sinks]::Vector{Sink}
    sinks = sinks[1:n_sinks]

    # Estimate the densities of each sink

    ## Select the bandwidth for the estimation
    ## Uses Silverman's rule of thumb
    sink1 = sinks[begin]
    grain1 = sink1[begin]
    inner_percentile = 95 # Filter outliers; ignore values outside the inner percentile
    alpha = 0.9 * scale_bandwidth # bandwidth "alpha" smooths the density estimate, 0.9 is the default

    bandwidth_matrix = zeros(length(sinks), length(grain1))
    for (i, sink) in enumerate(sinks)
        bandwidth_matrix[i, :] = default_bandwidth.(collect(eachmeasurement(sink)), alpha, inner_percentile)
    end

    bandwidths = dropdims(median(bandwidth_matrix, dims=1), dims=1) # column-wise median

    ## Note getmeasurements() gets the *names* of each measurement,
    ## whereas eachmeasurement() is an iterator for the *values* for each measurement

    ## Obtain the raw densities estimates
    ## The same measurement could (and likely!) have different supports for different sinks...
    raw_densities = make_densities.(sinks; bandwidths, inner_percentile)

    ## ...so we standardize them by resampling the densities on the same domain,
    ## which is the union of all intervels.
    densities, domains = standardize_KDEs(raw_densities; n_samples=discretization_size) #Made it to here so far

    ## We now have a list of densities and domain xs they are sampled at.
    ## Here, densities[i] are the densities for the ith sink (sinks[i]),
    ## and xs[j] is the domain for the jth measurement (densities[i][j]).
    ## Importantly, this is the same domain regardless of the sink!

    # Package the densities into a single order-3 tensor
    densitytensor = DensityTensor(densities, domains, sinks);
    setsourcename!(densitytensor, "sink");

    # Visualize the data in the tensor by plotting the densities for the first measurement
    p = plot_densities(densitytensor, "Age");
    plot!(xlabel="age (millions of years)",
        ylabel="probability density",
        legendtitle = "Sink",
        legend_columns=2,
        )
    display(p)

    # Find the best rank and perform the nonnegative decomposition Y=CF
    Y = copy(array(densitytensor)); # plain Array{T, 3} type for faster factorization

    Y_lateral_slices = eachslice(Y, dims=2)
    Y_lateral_slices .*= getstepsizes(densitytensor)

    Y = Y[:,1:n_features,:] # Only look at the first n_features features

    maxiter = 8000
    tol = 1e-5

    ranks = 1:size(Y)[1]
    Cs, Fs, all_rel_errors, final_errors, norm_grads, dist_Ncones = ([] for _ in 1:6)

    println("rank | n_iterations | final loss")
    for rank in ranks
        C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank; projection=:nnscale, maxiter, tol, rescale_Y=false);
        final_error = norm(Y - C*F) #rel_errors[end]
        push!.(
            (Cs, Fs, all_rel_errors, final_errors, norm_grads, dist_Ncones),
            (C, F, rel_errors, final_error, norm_grad, dist_Ncone)
        )
        @printf("%4i | %12i | %3.3g\n",
            rank, length(rel_errors), final_error)
    end
    final_rel_errors = convert(Vector{Float64}, final_errors)

    return final_rel_errors
end

# Parameters to check
scale_bandwidths = [0.5, 1, 1.5] # bandwidth parameter, larger -> more smoothing, 1 is the default (alpha=0.9)
proportion_of_grains_used = [1, 0.66, 0.33] # 75 grains total, so this uses 75, 50, 25 grains
# We use do it this way to make sure we take grains from every true source since the grains are ordered by source in the sinks

linecolors = [:orange, :green, :purple] # for each alpha
linestyles = [:solid, :dash, :dashdotdot] # for each proportion
lineoptions = [(:color => c, :linestyle => l) for c in linecolors, l in linestyles]

test_pairs = [(s, P) for s in scale_bandwidths, P in proportion_of_grains_used]

# Main computation
#final_errors = compute_final_loss_curve.(test_pairs)
final_errors_10 = map(f->f[1:10], final_errors)

# Final loss curve for all 9 tests
p = plot()
options = (:xlabel => "rank", legend_columns=3, linestyle=:dashdot, :color=>:blue, :linewidth=>2.5)
for (f, (s, P), extra) in zip(final_errors_10, test_pairs, lineoptions)
    plot!(f; ylabel="final loss", label="$((s, P))", options..., extra...)
end
display(p)

# Curvature for all 9 tests
p = plot(;leg=:bottomleft)
order = 5 # order to do the curature estimation
for (f, (s, P), extra) in zip(final_errors_10, test_pairs, lineoptions)
    plot!(standard_curvature(f; order); ylabel="standard curvature", label="$((s, P))", options..., extra...)
end
display(p)

# Show that we get rank = 3 each time
@show argmax.(standard_curvature.(final_errors_10; order))
@show all(x -> x == 3, argmax.(standard_curvature.(final_errors_10; order)) )

###############

# What happens if we use fewer features?

final_errors = compute_final_loss_curve.(1,1,1:7)
final_errors_10 = map(f->f[1:10], final_errors)

p = plot(; legend_title="# features")
labels = string.(1:7)
options = (:xlabel => "rank", legend_columns=3)#, linestyle=:dashdot, :color=>:blue, :linewidth=>2.5)
for (f, J) in zip(final_errors_10, labels)
    plot!(f; ylabel="final loss", label=J, options...)
end
display(p)

# Curvature for all 9 tests
p = plot(;leg=:bottomleft, legend_title="# features")
order = 4 # order to do the curature estimation
for (f, J) in zip(final_errors_10, labels)
    plot!(standard_curvature(f; order); ylabel="standard curvature",label=J, options...)
end
display(p)

@show argmax.(standard_curvature.(final_errors_10; order))

###############

# What happens if we use fewer sinks?

final_errors = compute_final_loss_curve.(1,1,7,3:20)
#final_errors_10 = map(f -> f[1:10], final_errors)

p = plot(; legend_title="# sinks")
labels = string.(3:20)
options = (:xlabel => "rank", legend_columns=3)#, linestyle=:dashdot, :color=>:blue, :linewidth=>2.5)
for (f, I) in zip(final_errors, labels)
    plot!(f; ylabel="final loss", label=I, options...)
end
display(p)

# Curvature for all 9 tests
p = plot(;leg=:bottomleft, legend_title="# sinks")
order = 3 # order to do the curature estimation
for (f, I) in zip(final_errors, labels)
    plot!(standard_curvature(f; order); ylabel="standard curvature",label=I, options...)
end
display(p)

@show argmax.(standard_curvature.(final_errors; order))

################


# What happens if we use fewer discretization sizes?
discretization_sizes = [16, 32, 64, 128]
final_errors = compute_final_loss_curve.(1,1,7,20,discretization_sizes)
final_errors_10 = map(f->f[1:10], final_errors)

p = plot(; legend_title="\$K\$")
labels = string.(discretization_sizes)
options = (:xlabel => "rank", legend_columns=3)#, linestyle=:dashdot, :color=>:blue, :linewidth=>2.5)
for (f, K) in zip(final_errors_10, labels)
    plot!(f; ylabel="final loss", label=K, options...)
end
display(p)

# Curvature for all 9 tests
p = plot(;leg=:bottomleft, legend_title="\$K\$")
order = 3 # order to do the curature estimation
for (f, K) in zip(final_errors_10, labels)
    plot!(standard_curvature(f; order); ylabel="standard curvature",label=K, options...)
end
display(p)
