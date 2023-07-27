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
C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank, tol=2, maxiter=1500; rescale_Y=false);

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
@show rel_error(C * F, Y) # < 1% error which is good!

# Repeat for rank in {1, 2, ..., 10}
# Note rank = 5 is the smallest rank where the relative error plateaus
ranks = 1:10
final_rel_error = zeros(length(ranks))
final_norm_grad = zeros(length(ranks))
final_dist_Ncone = zeros(length(ranks))
for rank in ranks
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank, tol=2, maxiter=1500; rescale_Y=false);
    final_rel_error[rank] = rel_errors[end]
    final_norm_grad[rank] = norm_grad[end]
    final_dist_Ncone[rank] = dist_Ncone[end]
    @show final_rel_error[rank], final_norm_grad[rank], final_dist_Ncone[rank]
end

# Visualizing convergence
## Can see that the relative error flatlines at rank = 5
options = (:label => false, :xlabel => "rank")
p = plot(final_rel_error; ylabel="final relative error", options...)
display(p)

## The optimal rank is the maximum curvature i.e. largest 2d derivative of the error
p = plot(d2_dx2(final_rel_error); ylabel="2nd derivative of relative error", options...)
display(p)

## Rank does not really have a corrolation with the final norm of gradient
p = plot(final_norm_grad; ylabel="final norm of gradient", options...)
display(p)

## Gradually increases because the parameter space is larger??
p = plot(final_dist_Ncone; ylabel="final distance to normal cone", options...)
display(p)
