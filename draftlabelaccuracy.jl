# Draft functions for checking the accuracy of labels of grains across all sinks/rocks

"""
    label_accuracy(labels, true_amounts::AbstractArray{T,3})

Calculates the number of correctly labeled grains and percentage of correct labels

Arguments
---------
- `labels`: Iterable of iterables with the labels for each grain (Vector{Vector{T}} or Tuple{Tuple})
- `true_amounts`: Number of grains from each source for every sink

Outputs
-------
- `n_correct_eachsink::Integer`
- `n_total_labels::Integer`
- `accuracy::Real``
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
