#= Short helpers for MTF.jl =#

"""Alias for an AbstractArray{T, 3}."""
const Abstract3Tensor{T} = AbstractArray{T, 3}

# TODO maybe use a custom symbol in the future
"""
    Base.*(A::AbstractMatrix, B::Abstract3Tensor)

Computes the Abstract3Tensor C where ``C_{ijk} = \\sum_{l=1}^L A_{il} * B_{ljk}``.

When the third dimention of B has length 1, this is equivilent to the usual
matrix-matrix multiplication. For this reason, we resuse the same symbol.

This is equivilent to the ``1``-mode product ``B \\times_1 A``.
"""
function Base.:*(A::AbstractMatrix, B::Abstract3Tensor)
    @einsum C[i,j,k] := A[i,l] * B[l,j,k]
    return C
end

"""
    combined_norm(u, v, ...)

Compute the combined norm of the arguments as if all arguments were part of one large array.

This is equivilent to `norm(cat(u, v, ...))`, but this
implimentation avoids creating an intermediate array.

```jldoctest
u = [3 0]
v = [0 4 0]
combined_norm(u, v)

# output

5.0
```
"""
combined_norm(vargs...) = sqrt(sum(_norm2, vargs))
_norm2(x) = norm(x)^2

"""
    ReLU(x)

Rectified linear unit; takes the max of 0 and x.
"""
ReLU(x) = max(0,x)

"""
    rel_error(x, xhat)

Compute the relative error between x (true value) and xhat (its approximation).

The relative error is given by:
```math
\\frac{\\lVert \\hat{x} - x \\rVert}{\\lVert x \\rVert}
```
See also [`mean_rel_error`](@ref).
"""
function rel_error(xhat, x)
    return norm(xhat - x) / norm(x)
end

"""
    mean_rel_error(X, Xhat; dims=(1,2))

Compute the mean relative error between the dims-order slices of X and Xhat.

The mean relative error is given by:
```math
\\frac{1}{N}\\sum_{j=1}^N\\frac{\\lVert \\hat{X}_j - X_j \\rVert}{\\lVert X_j \\rVert}
```
See also [`rel_error`](@ref).
"""
function mean_rel_error(Xhat, X; dims=(1,2))
    hatslices = eachslice(Xhat; dims)
    slices = eachslice(X; dims)
    return mean(@. norm(hatslices - slices) / (norm(slices)))
end

"""
    d2_dx2(y::AbstractVector{<:Real})

Approximates the 2nd derivative of a function using only given samples y of that function.

Assumes y came from f(x) where x was an evenly sampled, unit intervel grid.
Note the approximation uses centered three point finite differences for the
next-to-end points, and foward/backward three point differences for the begining/end points
respectively. The remaining interior points use five point differences.

Force only the third order approximation with third_order=true.
See [`d_dx`](@ref).
"""
function d2_dx2(y::AbstractVector{<:Real}, third_order::Bool=false)
    length(y) < 5 ? third_order=true : nothing
    return third_order ? _d2_dx2_3(y) : _d2_dx2_5(y)
end
# TODO is there a package that does this? The ones I've seen require the forward function.

function _d2_dx2_3(y::AbstractVector{<:Real})
    d = zero(y)
    for i in eachindex(y)[begin+1:end-1]
        d[i] = y[i-1] - 2*y[i] + y[i+1]
    end
    # Assume the same curvature at the end points
    d[begin] = d[begin+1]
    d[end] = d[end-1]
    return d
end

function _d2_dx2_5(y::AbstractVector{<:Real})
    d = zero(y)
    each_i = eachindex(y)

    # Interior Ppints
    for i in each_i[begin+2:end-2]
        d[i] = (-y[i-2] + 16*y[i-1] - 30*y[i] + 16*y[i+1] - y[i+2])/12
    end

    # Boundary and next-to boundary points
    i = each_i[begin+2]
    d[i-2] = (35*y[i-2] - 104*y[i-1] + 114*y[i] - 56*y[i+1] + 11*y[i+2])/12
    d[i-1] = (11*y[i-2] - 20*y[i-1] + 6*y[i] + 4*y[i+1] - y[i+2])/12

    i = each_i[end-2]
    d[i+1] = (-y[i-2] + 4*y[i-1] + 6*y[i] - 20*y[i+1] + 11*y[i+2])/12
    d[i+2] = (11*y[i-2] - 56*y[i-1] + 114*y[i] - 104*y[i+1] + 35*y[i+2])/12

    return d
end

"""
    d_dx(y::AbstractVector{<:Real})

Approximates the 1nd derivative of a function using only given samples y of that function.

Assumes y came from f(x) where x was an evenly sampled, unit intervel grid.
Note the approximation uses centered three point finite differences for the
next-to-end points, and foward/backward three point differences for the begining/end points
respectively. The remaining interior points use five point differences.

Force only the third order approximation with third_order=true.
See [`d2_dx2`](@ref).
"""
function d_dx(y::AbstractVector{<:Real}, third_order::Bool=false)
    length(y) < 5 ? third_order=true : nothing
    return third_order ? _d_dx_3(y) : _d_dx_5(y)
end

function _d_dx_3(y::AbstractVector{<:Real})
    d = zero(y)
    each_i = eachindex(y)

    for i in each_i[begin+1:end-1]
        d[i] = -y[i-1] + y[i+1]
    end

    i = each_i[begin+1]
    d[begin] = -3*y[i-1] + 4*y[i] - y[i+1]

    i = each_i[end-1]
    d[end] = y[i-1] - 4*y[i] + 3*y[i+1]
    return d
end

function _d_dx_5(y::AbstractVector{<:Real})
    d = zero(y)
    each_i = eachindex(y)

    # Interior Ppints
    for i in each_i[begin+2:end-2]
        d[i] = (2*y[i-2] - 16*y[i-1] + 16*y[i+1] - 2*y[i+2])/24
    end

    # Boundary and next-to boundary points
    i = each_i[begin+2]
    d[i-2] = (-50*y[i-2] + 96*y[i-1] - 72*y[i] + 32*y[i+1] - 6*y[i+2])/24
    d[i-1] = (-6*y[i-2] - 20*y[i-1] + 36*y[i] - 12*y[i+1] + 2*y[i+2])/24

    i = each_i[end-2]
    d[i+1] = (-2*y[i-2] + 12*y[i-1] - 36*y[i] + 20*y[i+1] + 6*y[i+2])/24
    d[i+2] = (6*y[i-2] - 32*y[i-1] + 72*y[i] - 96*y[i+1] + 50*y[i+2])/24

    return d
end

"""
    curvature(y::AbstractVector{<:Real})

Approximates the signed curvature of a function given evenly spaced samples.

Uses [`d_dx`](@ref) and [`d2_dx2`](@ref) to approximate the first two derivatives.
"""
function curvature(y::AbstractVector{<:Real})
    dy_dx = d_dx(y)
    dy2_dx2 = d2_dx2(y)
    return @. dy2_dx2 / (1 + dy_dx^2)^1.5
end

"""
    standard_curvature(y::AbstractVector{<:Real})

Approximates the signed curvature of a function, scaled to the unit box ``[0,1]^2``.

See [`curvature`](@ref).
"""
function standard_curvature(y::AbstractVector{<:Real})
    c = length(y) / maximum(y)
    dy_dx = c*d_dx(y)
    dy2_dx2 = c^2*d2_dx2(y)
    return @. dy2_dx2 / (1 + dy_dx^2)^1.5
end

"""

    projsplx(y::AbstractVector{<:Real})

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
