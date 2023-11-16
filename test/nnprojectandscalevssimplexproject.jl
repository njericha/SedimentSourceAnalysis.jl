# This is an example factorization that compares to tweaks to BCD subject to a simplex constraint
# 1) gradient step -> nonnegative project -> rescale
# 2) gradient step -> simplex project

using Random
using LinearAlgebra

"""
    BCD(Y, r, enforce_constraint!; maxiter=maxiter, tol=tol, seed=123)

Block coordinate descent using a particular function "enforce_constraint!"
to ensure B is in the simplex.

This solves the following problem:

    min ‖ Y - AB ‖  s.t.  A ≥ 0, B ≥ 0, ∑ Bᵢⱼ = 1.

Note we assume the correct rank "r" is given.

updateAB! is a function that says how to update A and B every iteration.
"""
function BCD(Y, r, updateAB!; maxiter=1000, tol=1e-2, seed=123)
    # Ensure we start from the same point if desired
    Random.seed!(seed)

    # Initialize
    m, n = size(Y)
    A = abs.(randn(m, r))
    B = abs.(randn(r, n))
    B ./= sum(B) # Since we want ∑ Bᵢⱼ = 1

    loss_pre_proj = zeros(maxiter)
    loss_post_proj = zeros(maxiter)

    i = 0

    # Ensure at least one iteration is performed
    while i == 0 || (norm(Y-A*B) > tol && i < maxiter)
        i += 1

        loss_pre_proj[i], loss_post_proj[i] = updateAB!(A, B, Y)
    end

    # Chop excess
    loss_pre_proj = loss_pre_proj[1:i]
    loss_post_proj = loss_post_proj[1:i]

    return A, B, loss_pre_proj, loss_post_proj
end

max_abs_eigval(A) = maximum(abs.(eigvals(A)))

function grad_stepA!(A, B, Y;w=1)
    BB = B*B'
    L = w*norm(BB)
    #L = 0.5*max_abs_eigval(BB)
    A .-= (A*BB .- Y*B') ./ L
end

function grad_stepB!(A, B, Y;w=1)
    AA = A'A
    L = w*norm(AA)
    #L = 0.5*max_abs_eigval(AA)
    B .-= (AA*B - A'Y) ./ L
end

ReLU(x) = max(0, x)
"""
    nnp!(A, B, Y)

Gradient step, then Nonnegative project, then jointly rescale
"""
function nnp!(A, B, Y;w=1)
    grad_stepA!(A, B, Y;w=w)
    A .= ReLU.(A)

    grad_stepB!(A, B, Y;w=w)
    loss_pre = norm(A*B - Y) # get loss before B is projected
    B .= ReLU.(B)

    # Jointly rescale
    B_sum = sum(B)
    B ./= B_sum
    A .*= B_sum

    loss_post = norm(A*B - Y)

    return loss_pre, loss_post
end

"""
    sxp!(A, B, Y)

Gradient step, then simplex project
"""
function sxp!(A, B, Y;w=1)
    grad_stepA!(A, B, Y;w=w)
    A .= ReLU.(A)

    grad_stepB!(A, B, Y;w=w)
    loss_pre = norm(A*B - Y) # get loss before B is projected
    B .= projsplx(B)

    loss_post = norm(A*B - Y)

    return loss_pre, loss_post
end

"""

    projsplx(y)

Projects (in Euclidian distance) the vector y into the simplex.

[1] Yunmei Chen and Xiaojing Ye, "Projection Onto A Simplex", 2011
"""
function projsplx(y)
    y_sorted = sort(y[:])
    n = length(y)
    i = n - 1
    t = 0 # need to ensure t has scope outside the while loop
    while true
        t = (sum(y_sorted[i+1:end]) - 1) / (n-i)
        if t >= y_sorted[i]
            break
        else
            i -= 1
        end

        if i >= 1
            continue
        else # i == 0
            t = (sum(y_sorted) - 1) / n
            break
        end
    end
    return ReLU.(y .- t)
end

"""
Same as projsplx, but optimized inner loop
"""
function projsplx2(y)
    y_sorted = sort(y[:])
    n = length(y)
    i = n - 1
    t = 0 # need to ensure t has scope outside the while loop
    while true
        t = (sum(y_sorted[i+1:end]) - 1) / (n-i)
        if t >= y_sorted[i]
            break
        elseif i == 1
            t = (sum(y_sorted) - 1) / n
            break
        end

        i -= 1
    end
    return ReLU.(y .- t)
end

###################
# Start of Script #
###################

Random.seed!(314) # need a different seed than what the algorithm will guess
m, r, n = 10, 5, 50

# Just consider matricies for now
A_true = abs.(randn(m, r))
B_true = abs.(randn(r, n))
B_true ./= sum(B_true) # Entries of B sum to one

Y_true = A_true*B_true
@assert size(Y_true) == (m, n)

# BCD algorithm Parameters
maxiter = 10000
tol = 1e-5

A_nnp, B_nnp, loss_pre_nnp, loss_post_nnp = BCD(Y_true, r, nnp!; maxiter, tol)
A_sxp, B_sxp, loss_pre_sxp, loss_post_sxp = BCD(Y_true, r, sxp!; maxiter, tol)

using Plots

p = plot(
    [loss_pre_nnp, loss_post_nnp, loss_pre_sxp, loss_post_sxp];
    labels = ["Pre NN Projection and Rescaling" "Post NN Projection and Rescaling" "Pre Simplex Projection" "Post Simplex Projection"],
    title = "2-norm Loss of Block Coord Descent",
    ylabel = "‖Y - AB‖",
    xlabel = "iteration",
    linewidth = 4,
    linestyle = [:solid :dash :solid :dash],
    yaxis=:log10)
display(p)

display(B_true)
println(sum(B_true),"\n")
display(B_nnp)
println(sum(B_nnp),"\n")
display(B_sxp)
println(sum(B_sxp),"\n")

@show norm(A_true*B_true - Y_true)
@show norm(A_nnp*B_nnp - Y_true)
@show norm(A_sxp*B_sxp - Y_true)

@show norm(B_true - B_nnp)
@show norm(B_true - B_sxp)
[5, 2, 1, 3, 4]
scatter(B_true[:], (B_nnp[[5, 2, 1, 3, 4],:])[:])
scatter(B_nnp[:], B_sxp[:])

@show norm(B_true - B_nnp)
@show norm(B_true - B_nnp[[5, 2, 1, 3, 4],:])


@show norm(B_true - B_sxp)
@show norm(B_true - B_sxp[[5, 2, 1, 3, 4],:])
