#=
example decomposition on a random tensor
=#
# set random seed for repeatability

# make random matrix and tensors A[i,r], B[r,j,k] for some R

# set Y[i,j,k] = A[i,r] * B[r,j,k]

# decompose Y

# compare learned A and B with the true A and B

include("../src/SedimentAnalysis.jl")
using .SedimentAnalysis
using Random
