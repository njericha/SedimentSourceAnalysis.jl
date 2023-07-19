"""
Import, store, and visualize sediment data.
"""
module SedimentTools
# TODO should import Base: convert, length etc. go here?
using LinearAlgebra
using DataFrames: dropmissing! # TODO see if this is better than the remove_missing function I made
using XLSX # TODO annotate exactly what is being used from each imported package
using OrderedCollections
using NamedArrays

#export ...

include("structs.jl")
include("importing.jl")
include("densityestimation.jl")
include("classifysource.jl")
include("viz.jl")
end
