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
"""
rel_error(xhat, x) = norm(xhat - x) / norm(x)
# TODO use Distances.jl to define more robust errors

"""
    d2_dx2(y::AbstractVector{<:Real})

Approximates the 2nd derivative of a function using only given samples y of that function.

Assumes y came from f(x) where x was an evenly sampled, unit intervel grid.
Note the approximation uses centered three point finite differences for the
next-to-end points, and foward/backward three point differences for the begining/end points
respectively. The remaining interior points use five point differences.

Force only the third order approximation with third_order=true.
"""
function d2_dx2(y::AbstractVector{<:Real}, third_order::Bool=false)
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
    for i in eachindex(y)[begin+2:end-2]
        d[i] = (-y[i-1] + 16*y[i-1] - 30*y[i] + 16*y[i+1] - y[i+2])/12
    end

    # Third order for the next-to-end points
    for i in eachindex(y)[[begin+1, end-1]]
        d[i] = y[i-1] - 2*y[i] + y[i+1]
    end

    # Assume the same curvature at the end points
    d[begin] = d[begin+1]
    d[end] = d[end-1]
    return d
end
