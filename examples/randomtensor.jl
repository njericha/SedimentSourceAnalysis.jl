#=
example decomposition on a random tensor
=#

using Random
using Statistics: mean
using MatrixTensorFactor
using SedimentSourceAnalysis
using Plots
using Printf
using Logging; disable_logging(Warn)

# Plot settings
Plots.resetfontsizes(); Plots.scalefontsizes(1.5)
plotfont="Computer Modern"
plotfontsize=13
default(legendfontsize=plotfontsize, plot_titlefontsize=plotfontsize+1, titlefont=plotfontsize, legendtitlefontsize=plotfontsize, fontfamily=plotfont, legendfont=plotfont)

# set random seed for repeatability
Random.seed!(314159265)

# make random matrix and tensors C[i,r], F[r,j,k] for some R
M, N, P = (50, 50, 50)
R = 5
C_true = abs.(randn((M, R)))
F_true = abs.(randn((R, N, P)))

# Normalize C. This gives us some scaling we can use for uniqueness
C_true ./= mean(C_true)

Y = C_true * F_true # set Y[i,j,k] = \sum_r C[i,r] * F[r,j,k]
# size(Y) == (M, N, P)

# Perform the nonnegative decomposition Y=CF
rank = 5
C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank, tol=3e-1, maxiter=2000; rescale_Y=false);

# Plot Convergence
plots = plot_convergence(rel_errors, norm_grad, dist_Ncone)
display.(plots)

# Normalize C. To match our known scaling. The relative sizes between entries do not change
mean_C = mean(C)
C ./= mean_C
F .*= mean_C

# Because the order of columns of C/ horizontal slices of F could be different,
# we permute them to match the orginial ordering. Now the solution (should be?) unique
match_sources!(C, F, C_true, F_true)

# Show the relative error between the learned and true matrix/tensor
# TODO use other metrics like RMSE and SNR
@show rel_error(C, C_true)
@show rel_error(F, F_true)
@show rel_error(C * F, Y) # < 1% error which is good!

# Repeat for rank in {1, 2, ..., 10}
ranks = 1:10
final_rel_error = zeros(length(ranks))
final_norm_grad = zeros(length(ranks))
final_dist_Ncone = zeros(length(ranks))
println("rank | n_iterations | relative error")
for rank in ranks
    C, F, rel_errors, norm_grad, dist_Ncone = nnmtf(Y, rank, tol=3e-1, maxiter=2000; rescale_Y=false);
    final_rel_error[rank] = rel_errors[end]
    final_norm_grad[rank] = norm_grad[end]
    final_dist_Ncone[rank] = dist_Ncone[end]
    @printf("%4i | %12i | %3.3g\n",
        rank, length(rel_errors), final_rel_error[rank])
end

# Visualizing convergence
## Can see that the relative error flatlines at rank = 5
options = (:label => false, :xlabel => "rank")
p = plot(final_rel_error; ylabel="final relative error", options...)
display(p)

## The optimal rank is the maximum curvature i.e. largest 2d derivative of the error
p = plot(standard_curvature(final_rel_error); ylabel="curvature of relative error", options...)
display(p)

## Rank does not really have a corrolation with the final norm of gradient
p = plot(final_norm_grad; ylabel="final norm of gradient", options...)
display(p)

## Gradually increases because the parameter space is larger??
p = plot(final_dist_Ncone; ylabel="final distance to normal cone", options...)
display(p)
