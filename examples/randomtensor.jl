#=
example decomposition on a random tensor
=#

using Random
using SedimentAnalysis

# set random seed for repeatability
Random.seed!(314159)

M, N, P = (50, 50, 50)
R = 5
# make random matrix and tensors C[i,r], F[r,j,k] for some R
C_true = abs.(randn((M, R)))
F_true = abs.(randn((R, N, P)))

Y = C_true * F_true # set Y[i,j,k] = C[i,r] * F[r,j,k]

# Perform the nonnegative decomposition Y=CF
rank = 5
C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank; rescale_Y=false);

# Plot Convergence
plots = plot_convergence(rel_errors, norm_grad, dist_Ncone)
display.(plots)

match_sources!(C, F, C_true, F_true) # Because the order of sources could be different

# Show the relative error between the learned and true matrix/tensor
# TODO use other metrics like RMSE and SNR
@show rel_error(C, C_true)
@show rel_error(F, F_true)
@show rel_error(C * F, Y)
