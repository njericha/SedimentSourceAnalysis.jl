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

using PyCall
using PlotlyJS: PlotlyJS, mgrid, isosurface, attr

using MatrixTensorFactor
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

# only take three measruements
selected_measurements = ["Eu_anomaly","Ti_temp", "Dy_Yb"]
#selected_measurements = ["Eu_anomaly","Th_U"]
#selected_measurements = ["Th_U","Ti_temp"]
sinks = [sink[selected_measurements] for sink in sinks]::Vector{Sink}
sink1 = sinks[1]

# Estimate the densities of each sink

## Select the bandwidth for the estimation
## Uses Silverman's rule of thumb
sink1 = sinks[begin]
inner_percentile = 95 # Filter outliers; ignore values outside the inner percentile
#alpha_ = 0.9 # bandwidth "alpha" smooths the density estimate, 0.9 is the default
             # this can denoise estimation

#bandwidth_matrix = zeros(length(sinks), length(grain1))
#for (i, sink) in enumerate(sinks)
#    bandwidth_matrix[i, :] = default_bandwidth.(collect(eachmeasurement(sink)), alpha_, inner_percentile)
#end

#bandwidths = dropdims(median(bandwidth_matrix, dims=1), dims=1) # column-wise median

#filter_inner_percentile()

# Make 3d KDEs
sink_to_measurements(sink) = [[grain[i] for grain in sink] for i in eachindex(selected_measurements)]

sink1_measurements = sink_to_measurements(sink1)

function inflate(a,b;amount=1.5)
    mean = 0.5*(a+b)
    return a - amount*(mean-a), b + amount*(b-mean)
end

fastkde = pyimport("fastkde")
M = 2^5 + 1
coords = [range(inflate(extrema(w)...)..., length = M) for w in sink1_measurements]

makekde(sink_measurements) = fastkde.pdf(sink_measurements...,axes = coords,do_xarray_subset=false, var_names = ["x", "y", "z"], num_points = M).values |> reversedims

PDF = fastkde.pdf(sink1_measurements...,axes = coords,do_xarray_subset=false, var_names = ["x", "y", "z"], num_points = M)
#axes = coords
# Extract PyObject to Julia arrays
reversedims(X) = permutedims(X, [i for i in length(size(X)):-1:1])
pdf = reversedims(PDF.values) #numpy uses row major ordering whereas julia uses column
# x_coord = PDF.x.values
# y_coord = PDF.y.values
# z_coord = PDF.z.values
x_coord, y_coord, z_coord = coords
XXX, YYY, ZZZ = mgrid(x_coord, y_coord, z_coord)
# plot 3d KDE
function plot3d(pdf3d; isomin=1e-5)
    PlotlyJS.plot(isosurface(
        x=XXX[:], # isosurface wants all entries in a single list
        y=YYY[:],
        z=ZZZ[:],
        value=pdf3d[:],
        opacity=0.6,
        isomin=isomin,
        isomax=maximum(pdf3d)*0.9,
        surface_count=5, # number of isosurfaces, 2 by default: only min and max
        colorbar_nticks=5, # colorbar ticks correspond to isosurface values
        caps=attr(x_show=false, y_show=false, z_show=false),
    )) |> display
end

plot3d(pdf;isomin=1e-5)

# Now do it for every sink

densities = [makekde(sink_to_measurements(sink)) for sink in sinks]
domains = [x_coord, y_coord, z_coord]

# Find the best rank and perform the nonnegative decomposition Y=CF
Y = cat(densities..., dims=4)
Y = permutedims(Y, [4,1,2,3])
#Y = copy(array(densitytensor)); # plain Array{T, 3} type for faster factorization
#grid_volume = eltype(Y)(x_coord.step * y_coord.step * z_coord.step) # this is Δx * Δy, area of the 2d grid

Y_slices = eachslice(Y, dims=1)
correction = sum.(Y_slices)
Y_slices ./= correction
#Y .*= grid_volume
#Y_lateral_slices = eachslice(Y, dims=2)
#Y_lateral_slices .*= getstepsizes(densitytensor)

@assert all(abs.(sum.(eachslice(Y, dims=1)) .- 1) .< 1e-5) # all horizontal slices should be within 1e-5 of 1.0

maxiter = 500

ranks = 1:4#size(Y)[1]
Cs, Fs, all_rel_errors, final_errors, norm_grads, dist_Ncones = ([] for _ in 1:6)

println("rank | n_iterations | final loss")
for rank in ranks
    #I, J, K, L = size(Y)
    #tol = 1e-5 / sqrt(rank*(I+J*K*L))
    tol = 1e-7
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank; projection=:nnscale, maxiter, tol, normalize=:slices, rescale_Y=false);
    final_error = norm(Y - C*F) #rel_errors[end]
    push!.(
        (Cs, Fs, all_rel_errors, final_errors, norm_grads, dist_Ncones),
        (C, F, rel_errors, final_error, norm_grad, dist_Ncone)
    )
    @printf("%4i | %12i | %3.3g\n",
        rank, length(rel_errors), final_error)
end
final_rel_errors = convert(Vector{Float64}, final_errors)
length.(all_rel_errors)
## The optimal rank is the maximum curvature i.e. largest 2d derivative of the error
options = (:label => false, :xlabel => "rank")
p = plot(final_rel_errors; ylabel="final loss", options...)
#plot(final_rel_errors; ylabel="relative error", linewidth=5, markershape=:circle, markersize=8, options...)
display(p)
order = 4
p = plot(d_dx(final_rel_errors;order); ylabel="derivative of final loss", options...)
display(p)
p = plot(d2_dx2(final_rel_errors;order); ylabel="2nd derivative of final loss", options...)
display(p)
p = plot(standard_curvature(final_rel_errors; order); ylabel="standard curvature\nof final loss", options...)
display(p)

## Extract the variables corresponding to the optimal rank
best_rank = 3 # argmax(standard_curvature(final_rel_errors))
@show best_rank
C, F, rel_errors, norm_grad, dist_Ncone = getindex.(
    (Cs, Fs, all_rel_errors, norm_grads, dist_Ncones),
    best_rank
)

# Rescale F to match the original scaling for densitytensor
F ./= grid_volume = eltype(Y)(x_coord.step * y_coord.step * z_coord.step)
factortensor = F
#F_lateral_slices = eachslice(F, dims=2)
#F_lateral_slices ./= getstepsizes(densitytensor)

#F_lateral_slices = eachslice(F_simplex, dims=2)
#F_lateral_slices ./= getstepsizes(densitytensor)

# Plot Convergence
plots = plot_convergence(rel_errors, norm_grad, dist_Ncone)
display.(plots)

# Compare learned C and F to the known sources
# Import data for ground truth F
filename = "./data/sundell2022/3Sources from Sundell et al 2022.xlsx"
sources = read_raw_data(filename)::Vector{Source}

## Confirm the measurements are the same and in the same order
sources = [source[selected_measurements] for source in sources]::Vector{Source}
@assert getmeasurements(sources[begin]) == getmeasurements(sinks[begin])

## Estimate densities for the known sources
## We pass in the previously calculated domains to ensure the densities are
## estimated on the same grid.
## They are wrapped in a tuple so the broadcasting is only on the sources.
# TODO make the make_densities2d function in densityestimation.jl, add domain and inner_percentile options

# true_densities = make_densities2d.(sources, (domains,); bandwidths, inner_percentile)

true_densities = [makekde(sink_to_measurements(source)) for source in sources]


#=
## Wrap in DensityTensor
factortensor_true = DensityTensor(true_densities, domains, getmeasurements(densitytensor));
setsourcename!(factortensor_true, "true source")
=#
factortensor_true = cat(true_densities..., dims=4)
factortensor_true = permutedims(factortensor_true, [4,1,2,3])
# Import data for ground truth C
source_filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
source_amounts = XLSX.readdata(source_filename,"Source proportions","B2:D21")
C_true = source_amounts / 75 # 75 Grains in each sink
coefficientmatrix_true = NamedArray(C_true, dimnames=("sink", "true source"))

# Ensure the factors in C and F are in the same order as the true densities
match_sources!(C, F, coefficientmatrix_true, factortensor_true)


for (r, density) in enumerate(eachslice(factortensor; dims=1))
    r == 3 ? nothing : continue
    plot3d(density) |> display
end

for (r, density) in enumerate(eachslice(factortensor_true; dims=1))
    r == 3 ? nothing : continue
    plot3d(density) |> display
end

indx = factortensor_true .> 1e-3
p = scatter(factortensor_true[indx], factortensor[indx];
    xlabel="true",
    ylabel="learned",
    label="(true, learned)",
    )
xy_line = collect(extrema(factortensor_true[:]))
plot!(xy_line, xy_line, label="y=x line")
display(p)
@show mean_rel_error(factortensor, factortensor_true; dims=1)
#@show mean_rel_error(C * factortensor, Y ./ grid_volume; dims=1)
#mean_rel_error(C * factortensor, Y; dims=1)
# Now go back and classify each grain as source 1, 2, or 3 based on the learned sources

## We should see a nice step pattern since the sinks have grains order by the source they
## came from. You would not expect to have this perfect ordering for real data.
## Idealy the "misses" have smaller likelihoods
## i.e. we are less confident which source the grain came from.

## For each grain, get the source index estimate, and the list of likelihoods for each source
## Start with just the first sink
#domains = getdomains(factortensor) # extract domains and stepsizes once to speed up code
#stepsizes = getstepsizes(factortensor)
source_labels, source_likelihoods = zip(
    map(g -> estimate_which_nd_source(g, factortensor; domains, all_likelihoods=true), sinks[1])...)

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
    map(g -> estimate_which_nd_source(g, factortensor_true; domains, all_likelihoods=true), sinks[1])...)

## Sort the likelihoods (does not mutate source_likelihoods)
## and find the log of the max/2nd highest likelihood
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
label_grains(factortensor) = [map(g -> estimate_which_nd_source(g, factortensor; domains), sink) for sink in sinks]
source_labels = label_grains(factortensor)
n_correct_eachsink, n_total_labels, accuracy = label_accuracy(source_labels, source_amounts)

## Then the distributions from the sources
true_source_labels = label_grains(factortensor_true)
true_source_n_correct_eachsink, _, true_source_accuracy = label_accuracy(true_source_labels, source_amounts)

@show accuracy
@show true_source_accuracy

# Compare learned lable proportions, to the proportion/coefficient matrix

## Created learned proportion matrix C_label_proportions
sinki = 1
C_label_proportions = similar(C)
for (i, learned_source_proportions) in enumerate(eachrow(C))
    learned_gain_label_proportions_i =
        [count(source_labels[i] .== r) for r in 1:length(learned_source_proportions)] / 75 #number of grains in each sink
    C_label_proportions[i, :] = learned_gain_label_proportions_i
end

## Plot comparison between learned proportions and labeled proportions
heatmap(C; clims = (0,1)) |> display
heatmap(C_label_proportions; clims = (0,1)) |> display
p = scatter(C[:], C_label_proportions[:];
    xlabel="learned proportion matrix entries",
    ylabel="learned grain label proportions",
    label="(matrix entries, grain proportions)",)
plot!([0,1], [0,1]; label="y=x line")
display(p)

## Plot comparison between true proportions and labeled proportions
p = scatter(C_true[:], C_label_proportions[:];
    xlabel="true proportions",
    ylabel="learned grain label proportions",
    label="(true, grain proportions)",)
plot!([0,1], [0,1]; label="y=x line")
display(p)

## Calculate MAE
mean_abs(learned, truth) = mean(abs.(learned - truth))
@show mean_abs(C, C_true)
@show mean_abs(C_label_proportions, C)
@show mean_abs(C_label_proportions, C_true)
