#=
Use data from /data/sundell2022
=#

using XLSX
using NamedArrays
using Plots
using Printf
using Random
using Statistics: mean
using SedimentAnalysis

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
display(densitytensor[:, "Age", 1:1]) #just the first 1 sample

## ..or by sink
println("densitytensor Sink 1 (horizontal slice)")
display(densitytensor[1, :, 1]) #just the first sample

# Visualize the data in the tensor by plotting the densities for the first measurement
measurement_names = getmeasurements(densitytensor) # ["Age", ...]
p = plot_densities(densitytensor, "Age");
display(p)

# Find the best rank and perform the nonnegative decomposition Y=CF
Y = array(densitytensor); # plain Array{T, 3} type for faster factorization
maxiter = 5000
tol = 1e-5
#rank = 3
#C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank);

ranks = 1:length(getmeasurements(grain1))
Cs, Fs, all_rel_errors, norm_grads, dist_Ncones = ([] for _ in 1:5)

println("rank | n_iterations | relative error")
for rank in ranks
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank; maxiter, tol);
    push!.(
        (Cs, Fs, all_rel_errors, norm_grads, dist_Ncones),
        (C, F, rel_errors, norm_grad, dist_Ncone)
    )
    @printf("%4i | %12i | %3.3g\n",
        rank, length(rel_errors), rel_errors[end])
end

## The optimal rank is the maximum curvature i.e. largest 2d derivative of the error
options = (:label => false, :xlabel => "rank")
p = plot((map(x -> x[end],all_rel_errors)); ylabel="relative error", options...)
display(p)

p = plot(d2_dx2(map(x -> x[end],all_rel_errors)); ylabel="2nd derivative of relative error", options...)
display(p)

## Extract the variables corresponding to the optimal rank
best_rank = argmax(d2_dx2(map(x -> x[end],all_rel_errors)))
@show best_rank
C, F, rel_errors, norm_grad, dist_Ncone = getindex.(
    (Cs, Fs, all_rel_errors, norm_grads, dist_Ncones),
    best_rank
)

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

## Package F into a DensityTensor with the same domain and measurements as Y
## Each collection of measurements is no longer a sink and is now a "factor"
factortensor = DensityTensor(F, domains, getmeasurements(densitytensor))
setsourcename!(factortensor, "learned source")

## Package C into a NamedArray to label the dimentions
coefficientmatrix = NamedArray(C, dimnames=("sink", "learned source"))

# Visualize and compare the coefficientmatrix and factortensor
# Note display should be called after all plots have been finalized
options = (:xlabel => "source", :ylabel => "sink")
p = heatmap(coefficientmatrix; title="Learned Coefficients", clims = (0,1), options...);
display(p)
p = heatmap(coefficientmatrix_true; title="True Coefficients", clims = (0,1), options...);
display(p)

diff_coefficientmatrix = NamedArray(abs.(coefficientmatrix_true - coefficientmatrix),
    dimnames=("sink", "source"))
p = heatmap(diff_coefficientmatrix; title="Absolute Difference Between Learned and True Coefficients",
    clims = (0,0.15), options...);
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

# Show the relative error between the learned and true matrix/tensor
# ex. 0.15 relative error can be thought of as 85% similarity
# TODO use other metrics like RMSE and SNR
@show rel_error(coefficientmatrix, coefficientmatrix_true) # The warning that names are different is O.K.
@show rel_error(factortensor, factortensor_true)
@show rel_error(coefficientmatrix * factortensor, densitytensor)

# Now go back and classify each grain as source 1, 2, or 3 based on the learned sources

## We should see a nice step pattern since the sinks have grains order by the source they
## came from. You would not expect to have this perfect ordering for real data.
## Idealy the "misses" have smaller likelihoods
## i.e. we are less confident which source the grain came from.

## For each grain, get the source index estimate, and the list of likelihoods for each source
## Start with just the first sink
source_indexes, source_likelihoods = zip(
    map(g -> estimate_which_source(g, factortensor, all_likelihoods=true), sinks[1])...)

## Sort the likelihoods, and find the log of the max/2nd highest likelihood
sort!.(source_likelihoods, rev=true) # descending order
loglikelihood_ratios = [log10(s_likelihoods[1] / (s_likelihoods[2] + eps())) for s_likelihoods in source_likelihoods]

p = plot_source_index(
    collect(source_indexes), loglikelihood_ratios;
    title="Grains' Estimated Source and Log Likelihood Ratio"
)
display(p)

# Compare against classification using the true densities
source_indexes, source_likelihoods = zip(
    map(g -> estimate_which_source(g, factortensor_true, all_likelihoods=true), sinks[1])...)

## Sort the likelihoods, and find the log of the max/2nd highest likelihood
sort!.(source_likelihoods, rev=true) # descending order
loglikelihood_ratios = [log10(s_likelihoods[1] / (s_likelihoods[2] + eps())) for s_likelihoods in source_likelihoods]

p = plot_source_index(
    collect(source_indexes), loglikelihood_ratios;
    title="Grains' Estimated True Source and Log Likelihood Ratio"
)
display(p)

## Can also see that some grains are mislabelled using the true densities,
## Where there are 3 common grains that are mislablled using either true or learned densities
