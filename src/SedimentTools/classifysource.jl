"""
    match_sources!(C, F, C_true, F_true)

Permute sources in C and F to match the ground truth C_true and F_true.

Similarity is checked by finding the source that minimizes norm(c - c_true) where c is the
column of C.

# Parameters
- `double_check::Bool`: When true, repeat for F and F_true with their horizontal slices, and assert the ordering is the same as before.
"""
function match_sources!(
    C::AbstractMatrix,
    F::AbstractArray{T, 3},
    C_true::AbstractMatrix,
    F_true::AbstractArray{T, 3};
    double_check=false,
    ) where T <: Real
    # make list to store which true source the columns of C match
    n_factors = size(C_true)[2]
    true_ordering = zeros(Integer, n_factors)

    # loop over every column of C and find the best matching column of C_true
    for (i, c_true) ∈ enumerate(eachcol(C_true))
        _, i_true = findmin(c -> norm(c - c_true), eachcol(C))
        true_ordering[i] = i_true
    end
    @assert allunique(true_ordering)

    if double_check #repeat for F
        true_ordering2 = zeros(Integer, n_factors)
        for (i, slice_true) ∈ enumerate(eachslice(F_true, dims=1))
            _, i_true = findmin(s -> norm(s - slice_true), eachslice(F, dims=1))
            true_ordering2[i] = i_true
        end
        @assert true_ordering == true_ordering2
    end

    # Swap columns of C and horizontal slices of F to the new ordering
    C .= @view C[:,true_ordering]
    F .= @view F[true_ordering,:,:]
    return true_ordering
end

"""
    _find_subinterval(domain, value)

Find which subinterval the data lies in for each measurement

```jdocstest
julia> test_scale = 1:10
1:10

julia> find_subinterval(test_scale, 4.8)
5
```
"""
function _find_subinterval(domain, value)
    val, k = findmin(x -> abs(x - value), domain)
    return k
end

"""
    _estimate_prob(source, domains, grain)

Estimates the probability of observing grain from the single factor
This assumes the distributions are normalized in the sence that
they carry equal weight e.i. the row_sums are 1.

Note the distributions before rescaling have a row_sum equal to 1/stepsize
with stepsize being the intervel width used in scale.
"""
function _estimate_prob(grain::Grain, source, domains, stepsizes)
    measurement_probabilities = zeros(n_measurements)
    for (j, (grain_value, density, domain, stepsize)) ∈ enumerate(zip(grain, eachrow(source), domains, stepsizes))
        subinterval_k = _find_subinterval(domain, grain_value)
        density_k = density[subinterval_k] # Get the density for observing a measurement on that intervel
        measurement_probabilities[j] = density_k * stepsize # estimate the probability by 0th order area
        # TODO estimate the probability with 1st order trapizoid
    end
    ## Multiply the probabilities to get the overall probability of
    ## observing the data in that volume
    return prod(measurement_probabilities)
end

"""
    (source_index, likelyhoods) = estimate_which_source(grain::Grain, F::DensityTensor)

Returns the likelyhood and source index of the mostly likely factor
the grain vector came from.
"""
function estimate_which_source(grain::Grain, F::DensityTensor)
    getmeasurements(grain) == getmeasurements(F) ||
        ArgumentError("Grain and F don't have matching measurements")
    sources = eachsource(F)
    likelyhoods = zeros(length(sources))
    domains = getdomains(F)
    stepsizes = stepsizes(F)
    for (i, source) ∈ enumerate(sources)
        prob = estimate_prob(grain, source, domains, stepsizes)
        likelyhoods[i] = prob
    end
    return findmax(likelyhoods), likelyhoods #find the one with the highest probability
end
