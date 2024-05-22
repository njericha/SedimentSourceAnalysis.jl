"""
    match_sources!(C, F, C_true, F_true)

Permute sources in C and F to match the ground truth C\\_true and F\\_true.

Similarity is checked by finding the source that minimizes norm(c - c\\_true) where c is the
column of C.

# Parameters
- `double_check::Bool`: When true, repeat for F and F_true with their horizontal slices, and assert the ordering is the same as before.
"""
function match_sources!(
    C::AbstractMatrix,
    F::AbstractArray{T, N},
    C_true::AbstractMatrix,
    F_true::AbstractArray{T, N};
    double_check=false,
    ) where {T <: Real, N}
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
    F .= @view F[true_ordering, ((:) for _ in 1:(N-1))...]
    return true_ordering
end

"""
    _find_subinterval(domain, value)

Find which subinterval the data lies in for each measurement.

Assumes the domain is in increasing order. Returns the first intervel (1, 2) and if the
value is less than every element in the domain, returns the last intervel
(length(domain)-1 , length(domain)).

```jdocstest
julia> test_scale = 1:10
1:10

julia> find_subinterval(test_scale, 4.8)
(4, 5)
```
"""
function _find_subinterval(domain, value)
    left = findlast(domain .<= value)
    if isnothing(left)
        left = 1
    elseif left == length(domain)
        left -= 1
    end
    return left, left + 1
end

"""
    _estimate_prob(source, domains, grain)

Estimates the probability of observing grain from the single factor
This assumes the distributions are normalized in the sense that
they carry equal weight i.e. the row_sums are 1.

Note the distributions before rescaling have a row_sum equal to 1/stepsize
with stepsize being the intervel width used in scale.
"""
function _estimate_prob(grain::Grain, source, domains, stepsizes)
    n_measurements = length(getmeasurements(grain))
    measurement_probabilities = zeros(n_measurements)
    for (j, (grain_value, density, domain, stepsize)) ∈ enumerate(zip(grain, eachrow(source), domains, stepsizes))
        left, right = _find_subinterval(domain, grain_value)
        density_k = (density[left] + density[right])/2 # Get the density for observing a measurement on that intervel
        measurement_probabilities[j] = density_k * stepsize # estimate the probability by 0th order area
    end
    ## Multiply the probabilities to get the overall probability of
    ## observing the data in that volume
    return prod(measurement_probabilities)
end

function _estimate_2d_prob(grain::Grain, source, (x_domain, y_domain))
    xleft, xright = _find_subinterval(x_domain, grain[1])
    yleft, yright = _find_subinterval(y_domain, grain[2])
    # Perform a 4 point average of the 2d density function
    prob = 0.25 * sum(source[idx...] for idx in [[xleft, yleft], [xright, yleft], [xleft, yright], [xright, yright]])
    #measurement_probabilities = zeros(n_measurements)
    #for (j, (grain_value, density, domain, stepsize)) ∈ enumerate(zip(grain, eachrow(source), domains, stepsizes))
    #    left, right = _find_subinterval(domain, grain_value)
    #    density_k = (density[left] + density[right])/2 # Get the density for observing a measurement on that intervel
    #    measurement_probabilities[j] = density_k * stepsize # estimate the probability by 0th order area
    #end
    ## Multiply the probabilities to get the overall probability of
    ## observing the data in that volume
    return prob::T where T <: Real #prod(measurement_probabilities)
end

function _estimate_nd_prob(grain::Grain, source, domains)
    #n = length(size(source)) # order of the source array
    endpoints = (_find_subinterval(domain, feature) for (domain, feature) in zip(domains, grain))
    # Perform a 2^n point average of the nd density function?
    # this is exponentialy many points for high dimentions...
    #prob = 2^(-n) * sum(source[idx...] for idx in combinations(enpoints))
    # two point "diagonal" estimate
    left_idxs, right_idxs = zip(endpoints...)
    prob = 0.5 * (source[left_idxs...] + source[right_idxs...])
    return prob::T where T <: Real
end

"""
    estimate_which_source(grain::Grain, F::DensityTensor; kwargs...)

Returns the likelihood and source index of the mostly likely factor
the grain vector came from.

# Returns
- (Default) `source_index::Integer`: The index of the most likely source
- (when `max_likelihoods==true`) `(maxlikelihood, source_index)`: The most likely source and its likelihood
- (when `all_likelihoods==true`) `likelihoods::Vector{Real}`: Likelihood grain came from each source
- (when both are true) `((maxlikelihood, source_index), likelihoods)`
"""
function estimate_which_source(
    grain::Grain,
    F::DensityTensor;
    max_likelihoods=false,
    all_likelihoods=false,
    domains = getdomains(F),
    stepsizes = getstepsizes(F)
 )
    getmeasurements(grain) == getmeasurements(F) ||
        ArgumentError("Grain and F don't have matching measurements")
    sources = eachsource(F)
    likelihoods = zeros(length(sources))

    for (i, source) ∈ enumerate(sources)
        prob = _estimate_prob(grain, source, domains, stepsizes)
        likelihoods[i] = prob
    end

    if max_likelihoods && all_likelihoods
        return findmax(likelihoods), likelihoods
    elseif max_likelihoods
        return findmax(likelihoods)
    elseif all_likelihoods
        return findmax(likelihoods)[2], likelihoods
    else
        return findmax(likelihoods)[2]
    end
end

"""
    estimate_which_2d_source(grain::Grain, F; kwargs...)

Returns the likelihood and source index of the mostly likely factor
the grain vector came from.

Works on densities F where each horizontal slice is a discretized 2D KDE.

See [`estimate_which_source`](@ref).

# Returns
- (Default) `source_index::Integer`: The index of the most likely source
- (when `max_likelihoods==true`) `(maxlikelihood, source_index)`: The most likely source and its likelihood
- (when `all_likelihoods==true`) `likelihoods::Vector{Real}`: Likelihood grain came from each source
- (when both are true) `((maxlikelihood, source_index), likelihoods)`
"""
function estimate_which_2d_source(
    grain::Grain,
    F;
    max_likelihoods=false,
    all_likelihoods=false,
    domains,
 )
    getmeasurements(grain) == 2 ||
        ArgumentError("Grain and F don't have matching measurements")
    sources = eachslice(F; dims=1) #eachsource(F)
    likelihoods = zeros(length(sources))

    for (i, source) ∈ enumerate(sources)
        prob = _estimate_2d_prob(grain, source, domains)
        likelihoods[i] = prob
    end

    if max_likelihoods && all_likelihoods
        return findmax(likelihoods), likelihoods
    elseif max_likelihoods
        return findmax(likelihoods)
    elseif all_likelihoods
        return findmax(likelihoods)[2], likelihoods
    else
        return findmax(likelihoods)[2]
    end
end

# TODO combine estimate_which_source functions together

"""
    estimate_which_nd_source(grain::Grain, F; kwargs...)

Returns the likelihood and source index of the mostly likely factor
the grain vector came from.

Works on densities F where each first order slice is a discretized n-Dimentional KDE.

See [`estimate_which_source`](@ref) and [`estimate_which_2d_source`](@ref).

# Returns
- (Default) `source_index::Integer`: The index of the most likely source
- (when `max_likelihoods==true`) `(maxlikelihood, source_index)`: The most likely source and its likelihood
- (when `all_likelihoods==true`) `likelihoods::Vector{Real}`: Likelihood grain came from each source
- (when both are true) `((maxlikelihood, source_index), likelihoods)`
"""
function estimate_which_nd_source(
    grain::Grain,
    F;
    max_likelihoods=false,
    all_likelihoods=false,
    domains,
 )

    sources = eachslice(F; dims=1) #eachsource(F)
    likelihoods = zeros(length(sources))

    for (i, source) ∈ enumerate(sources)
        prob = _estimate_nd_prob(grain, source, domains)
        likelihoods[i] = prob
    end

    if max_likelihoods && all_likelihoods
        return findmax(likelihoods), likelihoods
    elseif max_likelihoods
        return findmax(likelihoods)
    elseif all_likelihoods
        return findmax(likelihoods)[2], likelihoods
    else
        return findmax(likelihoods)[2]
    end
end

"""
    label_accuracy(labels, true_amounts::AbstractArray{T,3})

Calculates the number of correctly labeled grains and percentage of correct labels

Arguments
---------
- `labels`: Iterable of iterables with the labels for each grain (`Vector{Vector{T}}` or `Tuple{Tuple}`)
- `true_amounts`: Number of grains from each source for every sink

Outputs
-------
- `n_correct_eachsink::Integer`
- `n_total_labels::Integer`
- `accuracy::Real`
"""
function label_accuracy(labels, true_amounts::AbstractMatrix{T}) where T <: Integer
    # Check valid arguments
    n_total_labels = sum(length.(labels))
    all(sum.(eachrow(true_amounts)) .== n_total_labels) ||
        ArgumentError("Number of labels does not match grains in true_amounts")

    true_labels = _labels_from.(eachrow(true_amounts))

    n_correct_eachsink = [count(l .== l_true) for (l, l_true) in zip(labels, true_labels)]
    n_total_correct_labels = sum(n_correct_eachsink)
    accuracy = n_total_correct_labels / n_total_labels * 100

    return n_correct_eachsink, n_total_labels, accuracy
end

function label_accuracy(labels, true_amounts::AbstractMatrix{Any})
    return label_accuracy(labels, convert(Matrix{Integer}, true_amounts))
end

"""
    _labels_from(integers::AbstractVector{Integer})

Given a list of integers, make a new list containing
integers[1] many 1s, integers[2] many 2s, etc.

Example
-------
_labels_from([3,4,2]) == [1 1 1 2 2 2 2 3 3]
"""
function _labels_from(integers::AbstractVector{T}) where T <: Integer
    N = sum(integers)
    v = fill(1, N)
    idx = integers[begin]
    for (i, n) ∈ enumerate(integers[begin+1:end]) # skip first since v is filled with 1s
        v[idx+1:idx+n] .= i+1
        idx += n
    end
    return v
end

"""
    confidence_score(source_likelihoods; sorted=false)

Scores each grain between [0, 1] on how confident it came from the most likely source.

The likelihoods should be sorted in desending order. If they are presorted, use the
keyword `sorted=true`.

Input
-----
`source_likelihoods`: A list-of-lists or NTuple of Vectors of the likelihood each grain came
from each source.
"""
function confidence_score(source_likelihoods; sorted=false)
    if sorted
        return [log10(p[1] / (max(p[2], eps(p[1])))) for p in source_likelihoods]
    else
        return confidence_score(sort.(source_likelihoods, rev=true); sorted=true)
    end
end
