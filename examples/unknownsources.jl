#=
Use data from /data/lee2021
=#

using XLSX
using NamedArrays
using Plots
using MatrixTensorFactor
using SedimentAnalysis
using Printf
using Random
using Logging; disable_logging(Warn)

# Plot settings
Plots.resetfontsizes(); Plots.scalefontsizes(1.5)
plotfont="Computer Modern"
plotfontsize=13
default(legendfontsize=plotfontsize, plot_titlefontsize=plotfontsize+1, titlefont=plotfontsize, legendtitlefontsize=plotfontsize, fontfamily=plotfont, legendfont=plotfont)

# set random seed for repeatability
# does not randomize data, but the initialization for nnmtf
Random.seed!(314159265)

# Import data from excel file
filename = "./data/lee2021/Lee et al 2021 All Measurements.xlsx"
#filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
sinks = read_raw_data(filename)::Vector{Sink}

## Look at a sink and grain
sink1 = sinks[begin]
grain1 = sink1[begin]
println("Sink 1")
display(sink1)
println("Grain 1 in Sink 1")
display(grain1)

## Here we see each measurement displayed as well as the sampled value
## To get the names of each measurement, use getmeasurements(.)
@show getmeasurements(grain1)

# Estimate the densities of each sink

## Select the bandwidth for the estimation
## Uses Silverman's rule of thumb
#n_density_samples = 2^7
sink1 = sinks[begin]
inner_percentile = 95 # Filter outliers; ignore values outside the inner percentile
alpha_ = 0.9 # bandwidth "alpha" smooths the density estimate, 0.9 is the default
             # this can denoise estimation
bandwidths = default_bandwidth.(collect(eachmeasurement(sink1)), alpha_, inner_percentile)

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
display(densitytensor[:, 1:1, 1:1]) #just the first sample

## ..or by sink
println("densitytensor Sink 1 (horizontal slice)")
display(densitytensor[1:1, :, 1:1]) #just the first sample

# Visualize the data in the tensor by plotting the densities for the first measurement
measurement_names = getmeasurements(densitytensor) # ["Ages", ...]
p = plot_densities(densitytensor, "Age");
display(p)

# Perform the nonnegative decomposition Y=CF
Y = copy(array(densitytensor)); # plain Array{T, 3} type for faster factorization

#Y_lateral_slices = eachslice(Y, dims=2)
#Y_lateral_slices .*= getstepsizes(densitytensor)
Y_fibres = eachslice(Y, dims=(1,2))
Y_fibres ./= sum.(Y_fibres)

ranks = 1:5#size(Y)[1]
maxiter = 6000
tol = 1e-5
Cs, Fs, all_rel_errors, norm_grads, dist_Ncones = ([] for _ in 1:5)

println("rank | n_iterations | relative error")
for rank in ranks
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank; projection=:nnscale, maxiter, tol, rescale_Y=false);
    push!.(
        (Cs, Fs, all_rel_errors, norm_grads, dist_Ncones),
        (C, F, rel_errors, norm_grad, dist_Ncone)
    )
    @printf("%4i | %12i | %3.3g\n",
        rank, length(rel_errors), rel_errors[end])
end

## The optimal rank is the maximum curvature i.e. largest 2d derivative of the error
options = (:label => false, :xlabel => "rank")
p = plot(standard_curvature(map(x -> x[end],all_rel_errors)); ylabel="curvature of relative error", options...)
display(p)

p = plot((map(x -> x[end],all_rel_errors)); ylabel="relative error", options...)
display(p)

## Extract the variables corresponding to the optimal rank
best_rank = 3 #argmax(standard_curvature(map(x -> x[end],all_rel_errors)))
@show best_rank
C, F, rel_errors, norm_grad, dist_Ncone = getindex.(
    (Cs, Fs, all_rel_errors, norm_grads, dist_Ncones),
    best_rank
)

# Plot Convergence
plots = plot_convergence(rel_errors, norm_grad, dist_Ncone)
display.(plots)

F = DensityTensor(F, domains, getmeasurements(densitytensor))
setsourcename!(F, "learned source")

## Package C into a NamedArray to label the dimentions
C = NamedArray(C, dimnames=("sink", "learned source"))

# Visualize and compare the coefficientmatrix and factortensor
# Note display should be called after all plots have been finalized
p = heatmap(C; title="Learned Coefficients");
display(p)
learned_source_plots = source_heatmaps(F; title="Learned Densities for ");
display.(learned_source_plots)

plots = measurement_heatmaps(F; title="Learned Densities for ");
display.(plots);

@show rel_error(C * F, Y)

# Now go back and classify each grain as source 1, 2, or 3 based on the learned sources

## We should see a nice step pattern since the sinks have grains order by the source they
## came from. You would not expect to have this perfect ordering for real data.
## Idealy the "misses" have smaller likelihoods
## i.e. we are less confident which source the grain came from.

## For each grain, get the source index estimate, and the list of likelihoods for each source

## Start with just the first sink
source_indexes, source_likelihoods = zip(
    map(g -> estimate_which_source(g, F, all_likelihoods=true), sinks[1])...)

## Sort the likelihoods, and find the log of the max/2nd highest likelihood
loglikelihood_ratios = confidence_score(source_likelihoods)

p = plot_source_index(
    collect(source_indexes), loglikelihood_ratios;
    title="Grains' Estimated Source and Log Likelihood Ratio"
)
display(p)

all(sum.(eachslice(F, dims=(1,2))) .≈ 1)
sum.(eachslice(C, dims=(1)))

source_identification_per_sink = []
    for (sink_number, sink) in zip(eachindex(sinks), sinks)
        source_indexes, source_likelihoods = zip(
            map(g -> estimate_which_source(g, F, all_likelihoods=true), sink)...)
        loglikelihood_ratios = confidence_score(source_likelihoods)
        source_identification = Dict("sources" => collect(source_indexes), "loglikelihood_ratios" => loglikelihood_ratios)
        push!(source_identification_per_sink, Dict("name" => "sink $sink_number", "data" => source_identification))
    end
