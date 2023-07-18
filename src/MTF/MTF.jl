"""
Matrix-Tensor Factorization
"""
module MTF
using LinearAlgebra
using Plots
using Random
using Einsum

export combined_norm, nnmtf, plot_factors, rel_error

include("utils.jl")
include("matrixtensorfactorize.jl")

end # MTF