"""
Import, store, and visualize sediment data.
"""
module SedimentTools
# Required packages
using LinearAlgebra: norm
using XLSX: readxlsx, sheetnames
using OrderedCollections: OrderedDict
using NamedArrays: NamedArray, NamedMatrix, NamedVector, setnames!, setdimnames!, dimnames
#using ReusePatterns: @forward, forward # need to add the regular function `forward` because it is called by the macro @forward
using KernelDensity: UnivariateKDE, kde, pdf
using Plots: plot, plot!, scatter
using Base: AbstractVecOrTuple, Generator

#using Pkg
#Pkg.add(url="https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl.git")
using MatrixTensorFactor: default_bandwidth, DEFAULT_ALPHA, filter_inner_percentile, filter_2d_inner_percentile

# Method extentions
using Base: getindex, names, show
using Plots: Plots, heatmap
using NamedArrays: NamedArrays#, NamedArray
using MatrixTensorFactor: MatrixTensorFactor, make_densities, make_densities2d

# Exports
export Grain, DensityTensor, Rock, Sink, Source # Types
export array, getdomain, getdomains, getsourcename, getsourcenames,  getmeasurements
export getstepsizes, namedarray, getsink, getsource # Getters
export normalize_density_sums!, normalize_density_sums, setsourcename! # Setters
export eachdensity, eachmeasurement, eachsink, eachsource # Iterators
include("structs.jl")

export read_raw_data
include("importing.jl")

#export make_densities # Method extention to Sink
include("densityestimation.jl")

export estimate_which_source, estimate_which_2d_source, label_accuracy, match_sources!, confidence_score
include("classifysource.jl")

export measurement_heatmaps, plot_densities, source_heatmaps, plot_convergence, plot_source_index
include("viz.jl")

end # module SedimentTools
