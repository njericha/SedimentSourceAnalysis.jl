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

# Can see the structure of the tensor by showing different slices

## This contains the first density sample for each sink and measurment
## The "density=(-103.994,...)" gives the "x" value each density was sampled at
## Most of these are zero since it is the left side of the distributions' support
println("densitytensor first density sample (frontal slice)")
display(densitytensor[:, :, 1:1])

## Similarly slice by measurement...
println("densitytensor Age density (lateral slice)")
display(densitytensor[:, "Age", 1:1]) #just the first 1 sample

## ..or by sink
println("densitytensor Sink 1 (horizontal slice)")
display(densitytensor[1, :, 1]) #just the first sample

# Visualize the data in the tensor by plotting the densities for the first measurement
measurement_names = getmeasurements(densitytensor) # ["Age", ...]
p = plot_densities(densitytensor, "Age");
plot!(xlabel="age (millions of years)",
    ylabel="probability density",
    legendtitle = "Sink",
    legend_columns=2,
    )
display(p)

# Vizualize the first sink's KDE for the first measurement (Age)
raw_data_ages_sink1 = [grain["Age"] for grain in sink1]
KDE_ages_sink1 = densitytensor[1, "Age", :] #sink 1
p = histogram(raw_data_ages_sink1;
    normalize=true,
    bins=range(0, 200, length=6),
    label="histogram",
    alpha=0.25,
    color=:black,
    )
plot!(getdomain(densitytensor, "Age"), KDE_ages_sink1;
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
    xlabel="age (millions of years)",
    ylabel="probability density",
    )
display(p)

# Vizualize the first sink's KDE and sample points
p = plot(getdomain(densitytensor, "Age"), KDE_ages_sink1;
label="continuous KDE",
color=:blue,
xlabel="age (millions of years)",
ylabel="probability density"
)
scatter!(getdomain(densitytensor, "Age"), KDE_ages_sink1;
marker=:circ,
markercolor=:blue,
label="KDE discretization")
display(p)

# Coarse Vizualize the first sink's KDE and sample points
p = plot(getdomain(densitytensor, "Age"), KDE_ages_sink1;
label="continuous KDE",
color=:blue,
xlabel="age (millions of years)",
ylabel="probability density"
)
scatter!(getdomain(densitytensor, "Age")[1:4:end], KDE_ages_sink1[1:4:end];
marker=:circ,
markercolor=:blue,
label="KDE discretization")
display(p)

# Coarse Vizualize the first sink's KDE and sample points
p = plot(getdomain(densitytensor, "Age"), KDE_ages_sink1;
label="continuous KDE",
color=:blue,
xlabel="age (millions of years)",
ylabel="probability density"
)

display(p)

n_boxes = 14
step = length(KDE_ages_sink1) ÷ n_boxes
y_coarse=KDE_ages_sink1[1:step:end]
x_coarse=getdomain(densitytensor, "Age")[1:step:end]

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

# Find the best rank and perform the nonnegative decomposition Y=CF
Y = copy(array(densitytensor)); # plain Array{T, 3} type for faster factorization

Y_lateral_slices = eachslice(Y, dims=2)
Y_lateral_slices .*= getstepsizes(densitytensor)

### Finding optimal rank with singular values of (unfolded) Y
n_sinks = size(Y)[1]
Y_mat = reshape(Y, n_sinks, :)
#σ = svdvals(Y_mat)
#plot(σ)
#partial_sum = [sum(σ[i:n_sinks].^2) .^ 0.5 for i in 2:n_sinks]
#push!(partial_sum, 0)
#plot(partial_sum)
###

# 1 |            5 | 0.777
#    2 |          337 | 0.506
#    3 |         3001 | 0.324
#    4 |         6591 | 0.282
#    5 |         7000 | 0.255
#    6 |         7000 | 0.229
#    7 |         7000 | 0.199
#    8 |         7000 | 0.179
#    9 |         7000 | 0.166
#   10 |         7000 | 0.148
#   11 |         7000 | 0.136
#   12 |         7000 | 0.123
#   13 |         7000 | 0.112
#   14 |         7000 | 0.0997
#   15 |         7000 | 0.0909
#   16 |         7000 | 0.0803
#   17 |         7000 | 0.0711
#   18 |         7000 | 0.0593
#   19 |         7000 | 0.0568
#   20 |         7000 | 0.0407

final_errors_max_i7000_tol_1e6 = [
    0.777,
    0.506,
    0.324,
    0.282,
    0.255,
    0.229,
    0.199,
    0.179,
    0.166,
    0.148,
    0.136,
    0.123,
    0.112,
    0.0997,
    0.0909,
    0.0803,
    0.0711,
    0.0593,
    0.0568,
    0.0407,
]
v = final_errors_max_i7000_tol_1e6


maxiter = 7000
tol = 1e-6#1e-5

ranks = 1:size(Y)[1]
Cs, Fs, all_rel_errors, final_errors, norm_grads, dist_Ncones = ([] for _ in 1:6)

println("rank | n_iterations | final loss")
for rank in ranks
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank; maxiter, tol, rescale_Y=false);
    final_error = norm(Y - C*F) #rel_errors[end]
    push!.(
        (Cs, Fs, all_rel_errors, final_errors, norm_grads, dist_Ncones),
        (C, F, rel_errors, final_error, norm_grad, dist_Ncone)
    )
    @printf("%4i | %12i | %3.3g\n",
        rank, length(rel_errors), final_error)
end
final_rel_errors = v#convert(Vector{Float64}, final_errors)
length.(all_rel_errors)
## The optimal rank is the maximum curvature i.e. largest 2d derivative of the error
options = (:label => false, :xlabel => "rank")
p = plot(final_rel_errors; ylabel="final loss", options...)
#plot(final_rel_errors; ylabel="relative error", linewidth=5, markershape=:circle, markersize=8, options...)
display(p)
order = 4
p = plot(d_dx(final_rel_errors;order); ylabel="derivative of final loss", options...)
display(p)
p = plot(d2_dx2(final_rel_errors;order=5); ylabel="2nd derivative of final loss", options...)
display(p)
p = plot(standard_curvature(final_rel_errors[1:20]; order); ylabel="standard curvature\nof final loss", options...)
plot!(standard_curvature(final_rel_errors[1:4]; order); label="4", options...)
plot!(standard_curvature(final_rel_errors[1:5]; order); label="5", options...)
plot!(standard_curvature(final_rel_errors[1:6]; order); label="6", options...)
plot!(standard_curvature(final_rel_errors[1:20]; order); label="7", options...)
display(p)
p = plot(final_rel_errors; ylabel="standard curvature\nof final loss", options...)


function smooth(y;α=1) # α=0 is no smoothing, α=1 is as three way average
    z = copy(y)
    b = 3 - 2*α
    for i in eachindex(z)[begin+1:end-1]
        z[i] = (α*y[i-1] + b*y[i] + α*y[i+1])/3
    end
    return z
end


function curvature2(y; Δx=1, kwargs...)
    dy_dx = d_dx(y; kwargs...) ./ Δx
    dy2_dx2 = d2_dx2(y; kwargs...) ./ Δx^2
    return @. dy2_dx2 / (1 + dy_dx^2)^1.5
end

function standard_curvature2(y; kwargs...)
    Δx = 1 / (length(y) - 1) # An interval 0:10 has length(0:10) = 11, but measure 10-0 = 10
    y_max = maximum(y)
    dy_dx = d_dx(y; kwargs...) / (Δx * y_max)
    dy2_dx2 = d2_dx2(y; kwargs...) / (Δx^2 * y_max)
    return @. dy2_dx2 / (1 + dy_dx^2)^1.5
end

p = plot(smooth(final_rel_errors); ylabel="standard curvature\nof final loss", options...)

display(p)
p = plot(d_dx(final_rel_errors |> smooth); ylabel="standard curvature\nof final loss", options...)
p = plot(d2_dx2(final_rel_errors .^ 2  |> smooth;order=3); ylabel="standard curvature\nof final loss", options...)

p = plot(standard_curvature2(final_rel_errors |> smooth;order=3); ylabel="derivative of final loss", options...)
display(p)

x = 0:10
Δx = x[2]-x[1]
y = @. x^3

y = y / 1000 # normalize

x = x / 10

plot(x,y)
plot(x, d_dx(y) ./ Δx)
plot(x, d2_dx2(y) ./ Δx^2)
plot(x, curvature2(y; Δx))

plot(x, standard_curvature(y))
plot(x, standard_curvature2(y))

## Extract the variables corresponding to the optimal rank
best_rank = 3 # argmax(standard_curvature(final_rel_errors))
@show best_rank
C, F, rel_errors, norm_grad, dist_Ncone = getindex.(
    (Cs, Fs, all_rel_errors, norm_grads, dist_Ncones),
    best_rank
)

C_simplex, F_simplex, rel_errors_simplex, norm_grad_simplex, dist_Ncone_simplex = nnmtf(Y, 3; projection=:simplex, maxiter, tol, rescale_Y=false)
Yhat = C*F
Yhat_simplex = C_simplex*F_simplex

p = plot(rel_errors_simplex[2:length(rel_errors)];
    yaxis=:log10,
    label="simplex",
    xlabel="iteration number",
    ylabel="mean relative error")
plot!(rel_errors[2:end], label="nonnegative & rescale")
display(p)

p = plot(dist_Ncone_simplex[2:length(rel_errors)];
    yaxis=:log10,
    label="simplex",
    xlabel="iteration number",
    ylabel="vector-set distance")
plot!(dist_Ncone[2:end], label="nonnegative & rescale")
display(p)

# Rescale F to match the original scaling for densitytensor
F_lateral_slices = eachslice(F, dims=2)
F_lateral_slices ./= getstepsizes(densitytensor)

F_lateral_slices = eachslice(F_simplex, dims=2)
F_lateral_slices ./= getstepsizes(densitytensor)

# Plot Convergence
plots = plot_convergence(rel_errors, norm_grad, dist_Ncone)
display.(plots)

# Compare learned C and F to the known sources
# Import data for ground truth F
filename = "./data/sundell2022/3Sources from Sundell et al 2022.xlsx"
sources = read_raw_data(filename)::Vector{Source}

## Confirm the measurements are the same and in the same order
@assert getmeasurements(sources[begin]) == getmeasurements(sinks[begin])

## Estimate densities for the known sources
## We pass in the previously calculated domains to ensure the densities are
## estimated on the same grid.
## They are wrapped in a tuple so the broadcasting is only on the sources.
true_densities = make_densities.(sources, (domains,); bandwidths, inner_percentile)
#true_densities = make_densities.(sources, (domains,); inner_percentile) # possibly different bandwidths

## Wrap in DensityTensor
factortensor_true = DensityTensor(true_densities, domains, getmeasurements(densitytensor));
setsourcename!(factortensor_true, "true source")

# Import data for ground truth C
source_filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
source_amounts = XLSX.readdata(source_filename,"Source proportions","B2:D21")
C_true = source_amounts / 75 # 75 Grains in each sink
coefficientmatrix_true = NamedArray(C_true, dimnames=("sink", "true source"))

# Ensure the factors in C and F are in the same order as the true densities
match_sources!(C, F, coefficientmatrix_true, factortensor_true)
match_sources!(C_simplex, F_simplex, coefficientmatrix_true, factortensor_true)

# Compare results between nnscale vs. simplex projection
p = scatter(C[:], C_simplex[:];
xlabel="nonnegative project & rescale",
        ylabel="simplex",
        label="(true, learned)")
        xy_line = collect(extrema(C[:]))
plot!(xy_line, xy_line, label="y=x line")
display(p)

p = scatter(F[:], F_simplex[:];
xlabel="nonnegative project & rescale",
    ylabel="simplex",
    label="(true, learned)")
    xy_line = collect(extrema(F[:]))
plot!(xy_line, xy_line, label="y=x line")
display(p)

p = scatter(Yhat[:], Yhat_simplex[:];
xlabel="nonnegative project & rescale",
    ylabel="simplex",
    label="(true, learned)")
    xy_line = collect(extrema(Yhat[:]))
plot!(xy_line, xy_line, label="y=x line")
display(p)

scatter(C[:], C_simplex[:]) |> display
scatter(F[:], F_simplex[:]) |> display
scatter(Yhat[:], Yhat_simplex[:]) |> display

## Package F into a DensityTensor with the same domain and measurements as Y
## Each collection of measurements is no longer a sink and is now a "factor"
factortensor = DensityTensor(F, domains, getmeasurements(densitytensor))
setsourcename!(factortensor, "learned source")

## Package C into a NamedArray to label the dimentions
coefficientmatrix = NamedArray(C, dimnames=("sink", "learned source"))

# Visualize and compare the coefficientmatrix and factortensor
# Note display should be called after all plots have been finalized
options = (:xlabel => "source", :ylabel => "sink")
p = heatmap(coefficientmatrix; clims = (0,1), options...);#title="Learned Proportions",
display(p)
p = heatmap(coefficientmatrix_true; clims = (0,1), options...); #title="True Proportions"
display(p)

diff_coefficientmatrix = NamedArray(abs.(coefficientmatrix_true - coefficientmatrix),
    dimnames=("sink", "source"))
p = heatmap(diff_coefficientmatrix; title="Absolute Difference Between Learned and True Coefficients",
    clims = (0,0.15), options...);
display(p)

mae = mean(diff_coefficientmatrix)

p = plot([0, 1], [0, 1], label="y=x line", ribbon=mae, color=palette(:default)[3],
fillalpha=0.3)
plot!(
[coefficientmatrix_true[:][1], coefficientmatrix_true[:][1]],[coefficientmatrix_true[:][1], coefficientmatrix[:][1]];
label="absolute error", color=palette(:default)[2])
for i in eachindex(coefficientmatrix_true[:])
    plot!(
[coefficientmatrix_true[:][i], coefficientmatrix_true[:][i]],[coefficientmatrix_true[:][i], coefficientmatrix[:][i]];
label=nothing, color=palette(:default)[2])
end
scatter!(coefficientmatrix_true[:], coefficientmatrix[:],
    xlabel="true proportions",
    ylabel="learned proportions",
    label="(true, learned)",
    color=palette(:default)[1])

display(p)

@show mean(diff_coefficientmatrix)

learned_source_plots = source_heatmaps(factortensor; title="Learned Densities for ");
true_source_plots = source_heatmaps(factortensor_true; title="True Densities for ");

# Alternate learned and true plot so you can compare similar plots side-by-side
for (l, t) in zip(learned_source_plots, true_source_plots)
    display(l)
    display(t)
end

plots = measurement_heatmaps(factortensor; title="Learned Densities for ");
display.(plots);

D = normalize_density_sums(factortensor_true)
max_density = maximum(D)
function all_source_heatmaps(D::DensityTensor; kwargs...)
    plots = []
    D = normalize_density_sums(D)
    measurements = getmeasurements(D)
    #domain_length = length(getdomains(D)[begin])
    for source in eachsource(D)
        h = heatmap(
            array(source);
            yticks=(eachindex(measurements), measurements),
            #xticks=([1, domain_length],["1", "$(size(array(source))[2])"]),
            xlabel="sample index k",
            yflip=true,
            clims = (0, max_density),
            kwargs...
            )
        push!(plots, h)
    end
    return plots
end
p_learned = all_source_heatmaps(factortensor)
p_true = all_source_heatmaps(factortensor_true)
for (l, t) in zip(p_learned, p_true)
    display(l)
    display(t)
end
# Plot learned vs true densities (in original units) by source
for (i,(true_source, learned_source)) in enumerate(zip(eachsource(factortensor_true),eachsource(factortensor)))
    p = scatter(true_source[:], learned_source[:];
        xlabel="true densities",
        ylabel="learned densities",
        label="(true, learned)",
        title= "Source $i")
    xy_line = collect(extrema(true_source[:]))
    plot!(xy_line, xy_line, label="y=x line")
    display(p)
end

# Plot learned vs true densities (normalized) by source
Δx = getstepsizes(densitytensor)
for (i,(true_source, learned_source)) in enumerate(zip(eachsource(factortensor_true),eachsource(factortensor)))
    true_source_normalized = true_source .* Δx
    learned_source_normalized = learned_source .* Δx
    p = scatter(true_source_normalized[:], learned_source_normalized[:];
        xlabel="true densities",
        ylabel="learned densities",
        label="(true, learned)",)
        #title= "Source $i")
    xy_line = collect(extrema(true_source_normalized[:]))
    plot!(xy_line, xy_line, label="y=x line")
    display(p)
end
true_densities_normalized = factortensor_true .* Δx'
learned_densities_normalized = factortensor .* Δx'
p = scatter(true_densities_normalized[:], learned_densities_normalized[:];
    xlabel="true densities",
    ylabel="learned densities",
    label="(true, learned)",)
    #title= "Source $i")
xy_line = collect(extrema(true_densities_normalized[:]))
plot!(xy_line, xy_line, label="y=x line")
display(p)

# Plot learned vs true densities by measurement
for (name,true_measurement, learned_measurement) in zip(measurement_names, eachmeasurement(factortensor_true),eachmeasurement(factortensor))
    p = scatter(true_measurement[:], learned_measurement[:];
        xlabel="true",
        ylabel="learned",
        label="(true, learned)",
        title= "Density of $name")
    xy_line = collect(extrema(true_measurement[:]))
    plot!(xy_line, xy_line, label="y=x line")
    display(p)
end

# Show the relative error between the learned and true matrix/tensor
# ex. 0.15 relative error can be thought of as 85% similarity
# TODO use other metrics like RMSE and SNR
@show rel_error(coefficientmatrix, coefficientmatrix_true) # The warning that names are different is O.K.

factortensor_normalized = factortensor .* Δx'
factortensor_true_normalized = factortensor_true .* Δx'
@show mean_rel_error(factortensor, factortensor_true)
@show mean_rel_error(factortensor_normalized, factortensor_true_normalized)
@show mean_rel_error(coefficientmatrix * factortensor, densitytensor)

# Now go back and classify each grain as source 1, 2, or 3 based on the learned sources

## We should see a nice step pattern since the sinks have grains order by the source they
## came from. You would not expect to have this perfect ordering for real data.
## Idealy the "misses" have smaller likelihoods
## i.e. we are less confident which source the grain came from.

## For each grain, get the source index estimate, and the list of likelihoods for each source
## Start with just the first sink
domains = getdomains(factortensor) # extract domains and stepsizes once to speed up code
stepsizes = getstepsizes(factortensor)
source_labels, source_likelihoods = zip(
    map(g -> estimate_which_source(g, factortensor; domains, stepsizes, all_likelihoods=true), sinks[1])...)

## Sort the likelihoods, and find the log of the max/2nd highest likelihood
sort!.(source_likelihoods, rev=true) # descending order
loglikelihood_ratios = [log10(s_likelihoods[1] / (s_likelihoods[2] + eps())) for s_likelihoods in source_likelihoods]

p = plot_source_index(
    collect(source_labels), loglikelihood_ratios;
    title="Grains' Estimated Source and Log Likelihood Ratio"
)
grain_index = 0
for n_grains in source_amounts[1,:][1:end-1]
    grain_index += n_grains
    plot!([grain_index+0.5, grain_index+0.5], [1, 3];color=:black,legend=false)
end
plot!(title="",
xlabel="grain index n",
ylabel="estimated source r",)
display(p)

# Compare against classification using the true densities
source_labels, source_likelihoods = zip(
    map(g -> estimate_which_source(g, factortensor_true; domains, stepsizes, all_likelihoods=true), sinks[1])...)

## Sort the likelihoods, and find the log of the max/2nd highest likelihood
loglikelihood_ratios = confidence_score(source_likelihoods)

p = plot_source_index(
    collect(source_labels), loglikelihood_ratios;
    title="Grains' Estimated True Source and Log Likelihood Ratio"
)
grain_index = 0
for n_grains in source_amounts[1,:][1:end-1]
    grain_index += n_grains
    plot!([grain_index+0.5, grain_index+0.5], [1, 3];color=:black,legend=false)
end
plot!(title="",
xlabel="grain index n",
ylabel="estimated source r",)
display(p)

## Can also see that some grains are mislabelled using the true densities,
## Where there are 3 common grains that are mislablled using either true or learned densities

# Label every grain in every sink
## Use the learned distributions first
label_grains(factortensor) = [map(g -> estimate_which_source(g, factortensor; domains, stepsizes), sink) for sink in sinks]
source_labels = label_grains(factortensor)
n_correct_eachsink, n_total_labels, accuracy = label_accuracy(source_labels, source_amounts)

## Then the distributions from the sources
true_source_labels = label_grains(factortensor_true)
true_source_n_correct_eachsink, _, true_source_accuracy = label_accuracy(true_source_labels, source_amounts)

@show accuracy
@show true_source_accuracy

## Can see similar accuracy which implies it is not the learned densities that are
## innaccurate, but about 10% of the grains are simply "unlikely" to be from the source they
## came from.
