#=
Use data from /data/sundell2022
=#

using SedimentAnalysis

# Import data from excel file
filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
sinks = read_raw_data(filename)::Vector{Sink}

# Estimate the densities of each sink

## Select the bandwidth for the estimation
## Uses Silverman's rule of thumb
sink1 = @view sinks[begin]
inner_percentile = 95 # Filter outliers; ignore values outside the inner percentile
alpha = 1.5 # smooth density estimate, 0.9 is the default
bandwidths = default_bandwidth.(eachmeasurment(sink1), alpha, inner_percentile)

## Obtain the raw densities estimates
## The same measurment could (and likely!) have different supports for different sinks...
raw_densities = make_densities.(sinks; bandwidths, inner_percentile)

## ...so we standardize them by resampling the densities on the same domain,
## which is the union of all intervels.
densities, domains = standardize_KDEs(raw_densities)

## We now have a list of densities and domain xs they are sampled at.
## Here, densities[i] are the densities for the ith sink (sinks[i]),
## and xs[j] is the domain for the jth measurement (densities[i][j]).
## Importantly, this is the same domain regardless of the sink!

# Package the densities into a single order-3 tensor
densitytensor = DensityTensor(densities, domains, sinks)
setsourcename!(factortensor, "sink")

# Visualize the tensor

# Perform the nonnegative decomposition Y=CF
Y = array(densitytensor) # plain Array{T, 3} type for faster factorization
rank = 3
C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank)

## Package F into a DensityTensor with the same domain and measurments as Y
## Each collection of measurments is no longer a sink and is now a "factor"
factortensor = DensityTensor(F; domain, measurments=measurments(densitytensor))
setsourcename!(factortensor, "learned source")

## Package C into a NamedMatrix to label the dimentions
coefficientmatrix = NamedMatrix(C, dimnames=("sink", "learned source"))

# Compare learned C and F to the known sources
# Import data for ground truth F
filename = "./data/sundell2022/3Sources from Sundell et al 2022.xlsx"
sources = read_raw_data(filename)::Vector{Source}

## Confirm the measurements are the same and in the same order
@assert measurments(sources[begin]) == measurments(sink[begin])

## Estimate densities for the known sources
## We pass in the previously calculated domains to ensure the densities are
## estimated on the same grid.
true_densities = make_densities.(sources; domains, bandwidths, inner_percentile)

## Wrap in DensityTensor
factortensor_true = DensityTensor(true_densities; domain, measurments=measurments(densitytensor))
setsourcename!(factortensor_true, "true source")

# Import data for ground truth C
filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
source_amounts = XLSX.readdata(source_filename,"Source proportions","B2:D21")
C_true = source_amounts / 75 # 75 Grains in each sink
coefficientmatrix_true = NamedMatrix(C_true, dimnames=("sink", "true source"))

# Ensure the factors in C and F are in the same order as the true densities
# Visualize coefficientmatrix and factortensor
