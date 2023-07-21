"""
Import, store, and visualize sediment data.
"""
module SedimentTools
# TODO should import Base: convert, length etc. go here?
using LinearAlgebra
#using DataFrames: dropmissing! # TODO see if this is better than the remove_missing function I made
using XLSX # TODO annotate exactly what is being used from each imported package
using OrderedCollections
using NamedArrays
using ReusePatterns
using KernelDensity

export Grain, DensityTensor, Rock, Sink, Source # Types
export array, domain, domains, measurments, nammedarray, sink, source # Getters
export setsourcename! # Setters
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
