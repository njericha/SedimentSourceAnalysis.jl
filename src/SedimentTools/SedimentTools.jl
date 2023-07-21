"""
Import, store, and visualize sediment data.
"""
module SedimentTools
# Required packages # TODO annotate exactly what is being used from each imported package
using LinearAlgebra
using Statistics: std, quantile
using XLSX: readxlsx, sheetnames
using OrderedCollections: OrderedDict
using NamedArrays: NamedArray, setnames!, setdimnames!
using ReusePatterns: ReusePatterns, @forward
using KernelDensity: UnivariateKDE, kde, pdf

# Method extentions
using Base: Base, getindex, names, (:*) #, convert
using Plots: Plots, heatmap

# Exports
export Grain, DensityTensor, Rock, Sink, Source # Types
export array, domain, domains, measurments, nammedarray, sink, source # Getters
export normalize_density_sums!, normalize_density_sums, setsourcename! # Setters
export eachdensity, eachmeasurment, eachsink, eachsource # Iterators
include("structs.jl")

export read_raw_data
include("importing.jl")

export DEFAULT_ALPHA, DEFAULT_N_SAMPLES # Constants
export default_bandwidth, make_densities, standardize_KDEs # Functions
include("densityestimation.jl")

include("classifysource.jl")

include("viz.jl")
end
