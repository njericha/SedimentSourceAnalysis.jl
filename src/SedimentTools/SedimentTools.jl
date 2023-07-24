"""
Import, store, and visualize sediment data.
"""
module SedimentTools
# Required packages
using Statistics: std, quantile
using XLSX: readxlsx, sheetnames
using OrderedCollections: OrderedDict
using NamedArrays: NamedArray, NamedMatrix, NamedVector, setnames!, setdimnames!
using ReusePatterns: @forward, forward # need to add the regular function `forward` because it is called by the macro @forward
using KernelDensity: UnivariateKDE, kde, pdf
using Plots: plot

# Method extentions
using Base: getindex, names #, convert
using Plots: Plots, heatmap
using NamedArrays: NamedArrays, NamedArray

# Exports
export Grain, DensityTensor, Rock, Sink, Source # Types
export array, domain, domains, getsourcename, measurements, nammedarray, sink, source # Getters
export normalize_density_sums!, normalize_density_sums, setsourcename! # Setters
export eachdensity, eachmeasurement, eachsink, eachsource # Iterators
include("structs.jl")

export read_raw_data
include("importing.jl")

export DEFAULT_ALPHA, DEFAULT_N_SAMPLES # Constants
export default_bandwidth, make_densities, standardize_KDEs # Functions
include("densityestimation.jl")

export estimate_which_source, match_sources!
include("classifysource.jl")

export measurement_heatmaps, plot_densities, source_heatmaps
include("viz.jl")

end # module SedimentTools
