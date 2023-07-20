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
raw_densities = make_densities.(sinks, bandwidths, inner_percentile)

## ...so we standardize them by resampling the densities on the same domain,
## which is the union of all intervels.
densities, domains = standardize_KDEs(raw_densities)

## We now have a list of densities and domain xs they are sampled at.
## Here, densities[i] are the densities for the ith sink (sinks[i]),
## and xs[j] is the domain for the jth measurement (densities[i][j]).
## Importantly, this is the same domain regardless of the sink!

# Package the densities into a single order-3 tensor
densitytensor = DensityTensor(densities, domains, sinks)
