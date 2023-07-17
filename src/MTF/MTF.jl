"""
Matrix-Tensor Factorization
"""
module MTF
using Random
using LinearAlgebra
using Einsum

export rel_error

include("utils.jl")
include("matrixtensorfactorize.jl")

end # MTF
