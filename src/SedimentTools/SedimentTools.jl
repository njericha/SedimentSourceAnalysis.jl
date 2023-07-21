"""
Import, store, and visualize sediment data.
"""
module SedimentTools
# Required packages
using LinearAlgebra
using Statistics: std, quantile
using XLSX: readxlsx, sheetnames
using OrderedCollections: OrderedDict
using NamedArrays: NamedArray, NamedMatrix, NamedVector, setnames!, setdimnames!
using ReusePatterns: @forward, forward # need to add the regular function `forward` because it is called by the macro @forward
using KernelDensity: UnivariateKDE, kde, pdf
using Plots: heatmap

# Method extentions
using Base: getindex, names #, convert

# Exports
export Grain, DensityTensor, Rock, Sink, Source # Types
export array, domain, domains, getsourcename, measurments, nammedarray, sink, source # Getters
export normalize_density_sums!, normalize_density_sums, setsourcename! # Setters
export eachdensity, eachmeasurment, eachsink, eachsource # Iterators
include("structs.jl")

export read_raw_data
include("importing.jl")

export DEFAULT_ALPHA, DEFAULT_N_SAMPLES # Constants
export default_bandwidth, make_densities, standardize_KDEs # Functions
include("densityestimation.jl")

include("classifysource.jl")

export measurment_heatmaps, source_heatmaps
include("viz.jl")
end # module SedimentTools
