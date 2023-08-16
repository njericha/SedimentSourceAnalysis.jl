"""
Matrix-Tensor Factorization
"""
module MTF
using LinearAlgebra: norm
using Plots: plot
using Statistics: mean, median
using Random: randn
using Einsum: @einsum

# Method extentions
using Base: *

export Abstract3Tensor # Types
export combined_norm, d2_dx2, dist_to_Ncone, nnmtf, plot_factors, rel_error, mean_rel_error # Functions

include("utils.jl")
include("matrixtensorfactorize.jl")

end # MTF
