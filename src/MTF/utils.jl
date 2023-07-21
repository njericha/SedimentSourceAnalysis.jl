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

This is equivilent to `norm(cat(u, v, ...))``, but this
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
\\frac{\\lVert x - \\hat{x} \\rVert}{\\lVert x \\rVert}
```
"""
rel_error(x, xhat) = norm(x - xhat) / norm(x)
# TODO use Distances.jl to define a more robust error
