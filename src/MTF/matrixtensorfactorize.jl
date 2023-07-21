#=
Holds the block coordinate decent factorization function
and related helpers
=#

"""
    nnmtf(Y::Abstract3Tensor, R::Integer; kwargs...)

Non-negatively matrix-tensor factorizes an order 3 tensor Y with a given "rank" R.

Factorizes ``Y \\approx C F`` where ``\\displaystyle Y[i,j,k] \\approx \\sum_{r=1}^R C[i,r]*F[r,j,k]``
and the factors ``C, F \\geq 0`` are nonnegative.

Note there may NOT be a unique optimal solution

# Arguments
- `Y::Abstract3Tensor`: tensor to factorize
- `R::Integer`: rank to factorize Y (size(C)[2] and size(F)[1])

# Keywords
- `maxiter::Integer=100`: maxmimum number of iterations
- `tol::Real=1e-3`: desiered tolerance for the -gradient's distance to the normal cone
- `rescale::Bool=false`: scale F at each iteration so that the 3-fiber sums are 1 (on average). This also preprocesses the input `Y` in a similar way, and rescales the final F so Y=CF.
- `plot_F::Integer=0`: if not 0, plot F every plot_F iterations
- `names::AbstractVector{String}=String[]`: names of the slices of F to use for ploting

# Returns
- `C::Matrix{Float64}`: the matrix C in the factorization Y ≈ C * F
- `F::Array{Float64, 3}`: the tensor F in the factorization Y ≈ C * F
- `rel_errors::Vector{Float64}`: relative errors at each iteration
- `norm_grad::Vector{Float64}`: norm of the full gradient at each iteration
- `dist_Ncone::Vector{Float64}`: distance of the -gradient to the normal cone at each iteration
"""
function nnmtf(
    Y::Abstract3Tensor,
    R::Integer;
    maxiter::Integer=100,
    tol::Real=1e-3,
    rescale::Bool=true,
    plot_F::Integer=0,
    names::AbstractVector{String}=String[],
)
    # Extract Dimentions
    M, N, P = size(Y)

    # Initialize C, F
    C = abs.(randn((M, R)))
    F = abs.(randn((R, N, P)))

    problem_size = R*(M + N*P)

    # Scale Y if desired
    if rescale
        Y_input = copy(Y)
        Y, avg_factor_sums = rescaleY(Y)
    end

    # Initialize Looping
    i = 1
    rel_errors = zeros(maxiter)
    norm_grad = zeros(maxiter)
    dist_Ncone = zeros(maxiter)

    # Calculate initial relative error and gradient
    rel_errors[i] = rel_error(Y, C*F)
    grad_C, grad_F = calc_gradient(C, F, Y)
    norm_grad[i] = combined_norm(grad_C, grad_F)
    dist_Ncone[i] = dist_to_Ncone(grad_C, grad_F, C, F)

    # Convergence criteria. We "normalize" the distance vector so the tolerance can be
    # picked independent of the dimentions of Y and rank R
    converged(dist_Ncone, i) = dist_Ncone[i]/sqrt(problem_size) < tol

    # Main Loop
    # Ensure at least 1 step is performed
    while (i == 1) || (!converged(dist_Ncone, i) && (i < maxiter))
        if (plot_F != 0) && ((i-1) % plot_F == 0)
            plot_factors(F, names, appendtitle=" at i=$i")
        end

        updateC!(C, F, Y)
        updateF!(C, F, Y)

        rescale ? rescaleCF!(C, F) : nothing

        # Calculate relative error and norm of gradient
        i += 1
        rel_errors[i] = rel_error(Y, C*F)
        grad_C, grad_F = calc_gradient(C, F, Y)
        norm_grad[i] = combined_norm(grad_C, grad_F)
        dist_Ncone[i] = dist_to_Ncone(grad_C, grad_F, C, F)
    end

    # Chop Excess
    keep_slice = 1:i
    rel_errors = rel_errors[keep_slice]
    norm_grad = norm_grad[keep_slice]
    dist_Ncone = dist_Ncone[keep_slice]

    # Rescale F back if Y was initialy scaled
    if rescale
        # Compare:
        # If F_rescaled := avg_factor_sums * F,
        # Y_input ≈ C * F_rescaled
        #       Y ≈ C * F (Here, Y and F have normalized fibers)
        F .= avg_factor_sums * F
    end

    return C, F, rel_errors, norm_grad, dist_Ncone
end

"""
    dist_to_Ncone(grad_C, grad_F, C, F)

Calculate the distance of the -gradient to the normal cone of the positive orthant.
"""
function dist_to_Ncone(grad_C, grad_F, C, F)
    grad_C_restricted = grad_C[(C .> 0) .|| (grad_C .< 0)]
    grad_F_restricted = grad_F[(F .> 0) .|| (grad_F .< 0)]
    return sqrt(norm(grad_C_restricted)^2 + norm(grad_F_restricted)^2)
end

# TODO move this ploting function to SedimentTools? Or seperate viz.jl file?
"""
    plot_factors(F, names; appendtitle="")

Plot each horizontal slice of F. Names give the name of each vertical slice.
"""
function plot_factors(F, names=string.(eachindex(F[1,:,1])); appendtitle="")
    size(F)[2] == length(names) || ArgumentError("names should have the same length as size(F)[2]")
    fiber_sums = sum.(eachslice(F,dims=(1,2)))
    avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
    F_temp = avg_factor_sums^(-1) * F
    for (j, F_slice) ∈ enumerate(eachslice(F_temp,dims=1))
        p = heatmap(F_slice,
            yticks=(eachindex(F_slice[:,1]), names),
            xticks=([1, length(F_slice[1,:])],["0", "1"]),
            xlabel="Normalized Range of Values",
            title = "Learned Distributions for Factor $j" * appendtitle,
        )
        display(p)
    end
end

function updateC!(C, F, Y)
    @einsum FF[i,j] := F[i,p,q]*F[j,p,q]
    @einsum GG[i,j] := Y[i,p,q]*F[j,p,q]
    L = norm(FF) # TODO Optimize calculation to take advantage of symmetry of FF
    C .-= (C*FF .- GG) ./ L # gradient
    C .= ReLU.(C) # project
end

function updateF!(C, F, Y)
    CC = C'C
    L = norm(CC) # TODO Optimize calculation to take advantage of symmetry of CC
    F .-= (CC*F .- C'*Y) ./ L # gradient
    F .= ReLU.(F) # project
end

function calc_gradient(C, F, Y)
    @einsum FF[i,j] := F[i,p,q]*F[j,p,q]
    @einsum GG[i,j] := Y[i,p,q]*F[j,p,q]
    CC = C'C
    grad_C = C*FF .- GG
    grad_F = CC*F .- C'*Y
    return grad_C, grad_F
end

function rescaleCF!(C, F)
    fiber_sums = sum.(eachslice(F,dims=(1,2)))
    avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
    F .= avg_factor_sums^(-1) * F
    C .= C * avg_factor_sums
end

function rescaleY(Y)
    fiber_sums = sum.(eachslice(Y,dims=(1,2)))
    avg_factor_sums = Diagonal(mean.(eachrow(fiber_sums)))
    return avg_factor_sums^(-1) * Y, avg_factor_sums
end
