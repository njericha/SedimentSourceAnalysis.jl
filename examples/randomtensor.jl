#=
example decomposition on a random tensor
=#

using Random
using Statistics: mean
using SedimentAnalysis
using Plots

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
C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank, tol=2; rescale_Y=false);

# Plot Convergence
plots = plot_convergence(rel_errors, norm_grad, dist_Ncone)
display.(plots)

# Normalize C. This gives us a consistant scaling for C and F, and so
# C should match the true C. The relative sizes between entries do not change
mean_C = mean(C) / mean(C_true)
C ./= mean_C
F .*= mean_C

match_sources!(C, F, C_true, F_true) # Because the order of sources could be different

# Show the relative error between the learned and true matrix/tensor
# TODO use other metrics like RMSE and SNR
@show rel_error(C, C_true)
@show rel_error(F, F_true)
@show rel_error(C * F, Y)

# TODO repeat for rank in 1:10 and show how 5 is "optimal" in some sense
ranks = 1:10
final_rel_error = zeros(length(ranks))
final_norm_grad = zeros(length(ranks))
final_dist_Ncone = zeros(length(ranks))
for rank in ranks
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank, tol=2; rescale_Y=false);
    final_rel_error[rank] = rel_errors[end]
    final_norm_grad[rank] = norm_grad[end]
    final_dist_Ncone[rank] = dist_Ncone[end]
    @show final_rel_error[rank], final_norm_grad[rank], final_dist_Ncone[rank]
end

plot(final_rel_error) # Can see that the relative error flatlines at rank = 5
plot(final_norm_grad) # No corrolation with the final norm of gradient
plot(final_dist_Ncone) # Gradually increases because ??
