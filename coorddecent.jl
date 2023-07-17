#=
Holds the block coordinate decent factorization function
and related helpers
=#
using Random
using LinearAlgebra
using Einsum

ReLU(x) = max(0,x)

rel_error(x, xhat) = abs(x - xhat) / x

import Base: * # TODO use better style, maybe custom symbol
function *(A::AbstractMatrix{T}, B::Array{T,3}) where T
    @einsum C[i,j,k] := A[i,l] * B[l,j,k]
    return C
end

"""
    dist_to_Ncone_of_Rplus(grad_C, grad_F, C, F)

Calculate the distance of the -gradient to the normal cone of the positive orthant
"""
function dist_to_Ncone_of_Rplus(grad_C, grad_F, C, F)
    grad_C_restricted = grad_C[(C .> 0) .|| (grad_C .< 0)]
    grad_F_restricted = grad_F[(F .> 0) .|| (grad_F .< 0)]
    return sqrt(norm(grad_C_restricted)^2 + norm(grad_F_restricted)^2)
end

"""
Factorizes Y ≈ C F where ``Y[i,j,k] \\approx \\sum_r^R C[i,r]*F[r,j,k]``
and the factors C, F ≥ 0 are nonnegative

Note there may NOT be a unique optimal solution

If plot_progress is different from 0, plot the progress of F every plot_progress iterations
"""
function coorddecent(Y, R; maxiter=100, tol=1e-3, normalize_each_update = false, plot_progress = 0, names= nothing)
    # Extract Dimentions
    M, N, P = size(Y)

    # Initialize C, F
    C = abs.(randn((M, R)))
    F = abs.(randn((R, N, P)))

    problem_size = R*(M + N*P)

    # Create Updates
    function update_C!(C, F)
        @einsum FF[i,j] := F[i,p,q]*F[j,p,q]
        @einsum GG[i,j] := Y[i,p,q]*F[j,p,q]
        L = norm(FF) # TODO Optimize calculation to take advantage of symmetry of FF
        C .-= (C*FF .- GG) ./ L # gradient
        C .= ReLU.(C) # project
    end

    function update_F!(C, F)
        CC = C'C
        L = norm(CC) # TODO Optimize calculation to take advantage of symmetry of FF
        F .-= (CC*F .- C'*Y) ./ L # gradient
        F .= ReLU.(F) # project
    end

    function calc_gradient(C, F)
        @einsum FF[i,j] := F[i,p,q]*F[j,p,q]
        @einsum GG[i,j] := Y[i,p,q]*F[j,p,q]
        CC = C'C
        grad_C = C*FF .- GG
        grad_F = CC*F .- C'*Y
        return grad_C, grad_F
    end

    norm_gradient(grad_C, grad_F) = sqrt(norm(grad_C)^2 + norm(grad_F)^2)

    function rescaleCF!(C, F)
        fiber_sums = sum.(eachslice(F,dims=(1,2)))
        avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
        F .= avg_factor_sums^(-1) * F # TODO make more accurate scaling
        C .= C * avg_factor_sums
    end

    function plot_factors(F, i, names)
        fiber_sums = sum.(eachslice(F,dims=(1,2)))
        avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
        F_temp = avg_factor_sums^(-1) * F
        for (j, F_slice) ∈ enumerate(eachslice(F_temp,dims=1))
            p = heatmap(F_slice,
                yticks=(eachindex(F_slice[:,1]), names),
                xticks=([1, length(F_slice[1,:])],["0", "1"]),
                xlabel="Normalized Range of Values",
                title = "Learned Distributions for Factor $j at i=$i"
            )
            display(p)
        end
    end

    # Initialize Looping
    i = 1
    error = zeros(maxiter)
    norm_grad = zeros(maxiter)
    dist_Ncone = zeros(maxiter)

    # Calculate initial relative error and gradient
    error[i] = norm(Y - C*F) / norm(Y)
    grad_C, grad_F = calc_gradient(C, F)
    norm_grad[i] = norm_gradient(grad_C, grad_F)
    dist_Ncone[i] = dist_to_Ncone_of_Rplus(grad_C, grad_F, C, F)

    #not_converged(error, i) = rel_error(error[i], error[i-1]) > tol
    #not_converged(error, i) = error[i] > tol
    #not_converged(norm_grad, i) = (norm_grad[i]/norm_grad[1]) > tol # previous one I used
    not_converged(dist_Ncone, i) = dist_Ncone[i]/sqrt(problem_size) > tol # "normalize" the vector by √n
    #not_converged(dist_Ncone, i) = (dist_Ncone[i]/dist_Ncone[1]) > tol

    # Main Loop
    # Ensure at least 1 step is performed
    while (i == 1) || (not_converged(dist_Ncone, i) && (i < maxiter))
        if (plot_progress != 0) && ((i-1) % plot_progress == 0)
            plot_factors(F, i, names)
        end
        update_C!(C, F)
        update_F!(C, F)

        normalize_each_update ? rescaleCF!(C, F) : nothing

        # Calculate error and norm of gradient
        i += 1
        error[i] = norm(Y - C*F) / norm(Y)
        grad_C, grad_F = calc_gradient(C, F)
        norm_grad[i] = norm_gradient(grad_C, grad_F)
        dist_Ncone[i] = dist_to_Ncone_of_Rplus(grad_C, grad_F, C, F)
    end
    # Chop Excess
    keep_slice = 1:i
    error = error[keep_slice]
    norm_grad = norm_grad[keep_slice]
    dist_Ncone = dist_Ncone[keep_slice]
    return C, F, error, norm_grad, dist_Ncone
end

function onehot_rows!(X)
    for (i, row) ∈ enumerate(eachrow(X))
        max_val = findmax(row)[1]
        X[i,:] = (row .≈ max_val) # ex. [0, 0, 1, 0] where row = [2, -1, 4, 1]
    end
end

"""
Factorizes Y ≈ C F where Y[i,j] ≈ \\sum_r^R C[i,r]*F[r,j]
and the factors C, F ≥ 0 are nonnegative
and the coefficients C[i,j] ∈ {0, 1}

Note there may NOT be a unique optimal solution
"""
function binarycoorddecent(Y, R; maxiter=100, tol=1e-3)
    # Extract Dimentions
    M, N = size(Y)

    # Initialize C, F
    C = abs.(randn((M, R)))
    F = abs.(randn((R, N)))

    # Create Updates
    function update_C!(C, F)
        FF = F*F'
        L = norm(FF) # TODO Optimize calculation to take advantage of symmetry of FF
        C .-= (C*FF .- Y*F') ./ L # gradient
        onehot_rows!(C) # project
    end

    function update_F!(C, F)
        CC = C'C
        L = norm(CC) # TODO Optimize calculation to take advantage of symmetry of CC
        F .-= (CC*F .- C'*Y) ./ L # gradient
        F .= ReLU.(F) # project
        factor_norms = Diagonal(norm.(eachrow(F)))
        C[:,:] = C * factor_norms
        F[:,:] = factor_norms^(-1) * F
    end

    function norm_gradient(C, F)
        FF = F*F'
        CC = C'C
        grad_C = C*FF .- Y*F'
        grad_F = CC*F .- C'*Y
        return sqrt(norm(grad_C)^2 + norm(grad_F)^2)
    end

    # Initialize Looping
    i = 1
    error = zeros(maxiter)
    norm_grad = zeros(maxiter)

    # Calculate initial relative error and gradient
    error[i] = norm(Y - C*F) / norm(Y)
    norm_grad[i] = norm_gradient(C, F)

    #not_converged(error, i) = rel_error(error[i], error[i-1]) > tol
    #not_converged(error, i) = error[i] > tol
    not_converged(norm_grad, i) = (norm_grad[i]/norm_grad[1]) > tol

    # Main Loop
    # Ensure at least 1 step is performed
    while (i == 1) || (not_converged(norm_grad, i) && (i < maxiter))
        update_F!(C, F)
        update_C!(C, F)

        # Calculate error and norm of gradient
        i += 1
        error[i] = norm(Y - C*F) / norm(Y)
        norm_grad[i] = norm_gradient(C, F)
    end
    # Chop Excess
    error = error[1:i]
    norm_grad = norm_grad[1:i]
    return C, F, error, norm_grad
end
