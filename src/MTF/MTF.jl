"""
Matrix-Tensor Factorization
"""
module MTF
using LinearAlgebra # TODO annotate which functions are used
using Plots
using Random
using Einsum

# Method extentions
using Base: *

export Abstract3Tensor # Types
export combined_norm, nnmtf, plot_factors, rel_error # Functions

include("utils.jl")
include("matrixtensorfactorize.jl")

end # MTF
