#########
# Grain #
#########

"""
    Grain(v::AbstractVector{T}, measurement_names::AbstractVector{String})

Struct to hold grain level data
"""
const Grain = NamedVector{T} where T <: Real

"""
    getmeasurements(g::Grain)
    getmeasurements(s::Sink)
    getmeasurements(D::DensityTensor)

Getter for the measurement names.
"""
getmeasurements(g::Grain) = names(g, 1) # names(g) from NamedArray returns Vector{Vector{T}}

# Very hand-wavey stuff so you can call Grain(...)
"""
    Grain(v::AbstractVector{T}, measurement_names::AbstractVector{String})

Main constructor for a Grain.
"""
function (::Type{S})(v::AbstractVector{T}, measurement_names::AbstractVector{String}) where {S<:Grain, T<:Real}
    return NamedArray(v, (measurement_names,), ("measurement",))::Grain
end

#################
# Sinks / Rocks #
#################
"""
    Sink(grain1, grain2, ...)
    Sink([grain1, grain2, ...])

Collects a list of Grains into a Rock/Sink.

Ensures all Grains have the same names and are in the same order.
Struct to hold sink level data
"""
const Sink = Vector{Grain}

function Base.show(io::IO, x::MIME"text/plain", S::Sink)
    S_named = hcat(S...)
    setdimnames!(S_named, "grain", 2)
    println("$(length(S))-element Vector{Grain} with data:")
    show(io::IO, x::MIME"text/plain", S_named)
end

"""Gets the names of measurements from a Sink"""
getmeasurements(s::Sink) = iszero(length(s)) ? String[] : getmeasurements(s[1])

"""Gets all values of the measurement `k` in the Sink."""
Base.getindex(s::Sink, k::String) = collect(g[k] for g ∈ s)

"""Iterator for a list of values of each measurement"""
eachmeasurement(s::Sink) = (s[m] for m in getmeasurements(s))

"""
    Sink(grain1, grain2, ...)
    Sink([grain1, grain2, ...])

Collects a list of Grains into a Rock/Sink.

Ensures all Grains have the same names and are in the same order.
"""
function (::Type{S})(vec_of_grains::AbstractVector{Grain}) where S <: Sink # each element is a grain
    @assert allequal(getmeasurements.(vec_of_grains))
    return collect(vec_of_grains)::Sink
end
(::Type{S})(vec_of_grains::AbstractVector{Grain}...) where S <: Sink = Sink(vec_of_grains)

# Define aliases so Rock or Source can be used in place of Sink when those terms make more
# sense in those contexts.
"""Alias for [`Sink`](@ref)"""
const Rock = Sink
"""Alias for [`Sink`](@ref)"""
const Source = Sink

#################
# DensityTensor #
#################

# Idealy this would be a plain NamedArray, where we store the domain as
# the names along the third axis. But we have a different domain for each
# lateral slice j :( This means we can wrap the named array in a new type
# since Julia can only subtype abstract types...

# Attempted to use the @forward macro from ReusePatterns, but it is having issues with
# 1) packages that NamedArray extends methods but SedimentAnalysis does not and
# 2) breaks when the scope changes (something to do with where @eval is run)

# SO we will use a plain NamedArray where we hide the domains in the 3rd axis names

"""
    DensityTensor(KDEs, domains, sinks)
    DensityTensor(array, domains, measurement_names; kw...)

An order 3 array to hold the density distributions for multiple sinks.

KDEs is a Vector{Vector{Vector{T}}} like type whereas array is an Array{T,3} like type.

Call [`setsourcename!`](@ref) to set the source name (name of first dimention).
"""
const DensityTensor{T} = NamedArray{T, 3} where T <: Real

"""
    DensityTensor(array, domains, measurement_names; kw...)

Low level wrapper for DensityTensor.

Call [`setsourcename!`](@ref) to set the source name (name of first dimention).
"""
function (::Type{D})(
    array::AbstractArray{T, 3},
    domains::AbstractVector{<:AbstractVector{T}},
    measurement_names::AbstractVector{String};
    kw...
    ) where {D <: DensityTensor, T <: Real}

    # Make names for each index
    source_names = collect(eachindex(array[:,begin,begin]))
    density_names = collect(zip(domains...))
    names = (source_names, measurement_names, density_names)

    # Wrap in a named array
    densitytensor = NamedArray(array; names, kw...)

    # Set axis/dimention names
    setdimnames!(densitytensor, "measurement", 2)
    setdimnames!(densitytensor, "density", 3)

    return densitytensor::DensityTensor
end

"""Turn a DensityTensor into a NamedArray type."""
namedarray(D::DensityTensor) = D::NamedArray

"""Turn a DensityTensor or NamedArray into a plain Array type."""
array(D::DensityTensor) = D.array
array(N::NamedArray) = N.array

"""
    DensityTensor(KDEs, domains, sinks)

High level constructor which checks there is the correct number of sinks and measurements
"""
function (::Type{D})(
    KDEs::AbstractVecOrTuple{AbstractVecOrTuple{AbstractVecOrTuple{T}}},#UnivariateKDE
    domains::AbstractVector{<: AbstractVector{U}},
    sinks::AbstractVector{Sink},
    ) where {D <: DensityTensor, T <: Real, U <: Real}
    # Argument Handeling
    allequal(getmeasurements.(sinks)) ||
        ArgumentError("All sinks must have the same measurements in the same order.")
    length(sinks) == length(KDEs) ||
        ArgumentError("Must be the same number of sinks as there are lists of KDEs.")
    measurement_names = getmeasurements(sinks[begin])

    return DensityTensor(KDEs, domains, measurement_names)
end

"""
    DensityTensor(KDEs, domains, measurement_names)

Main high level constructor for DensityTensor
"""
function (::Type{D})(
    KDEs::AbstractVecOrTuple{AbstractVecOrTuple{AbstractVecOrTuple{T}}},#UnivariateKDE
    domains::AbstractVector{<: AbstractVector{U}},
    measurement_names::AbstractVector{String},
    ) where {D <: DensityTensor, T <: Real, U <: Real}
    # Argument Handeling
    length(measurement_names) == length(KDEs[begin]) ||
        ArgumentError("Must be the same number of measurements as there are KDEs for each sink.")
    allequal(length.(domains)) ||
        ArgumentError("Domains have different number of samples.")
    length(KDEs[begin][begin]) == length(domains[begin]) ||
        ArgumentError("Number of density samples in KDEs does not match number of domain samples.")

    # Magic line to turn the KDEs into an order-3 tensor
    # First two `cat`s take the list of list of list of T into an Array{T, 3} type
    # permutedims transposes the tensor to the right ordering of dimentions
    data = permutedims(cat(cat.(KDEs..., dims=2)..., dims=3), [2,3,1])

    # Confirm all the dimentions are in the right order
    @assert size(data) == (length(KDEs), length(measurement_names), length(domains[begin]))

    # Wrap again in a DensityTensor to store the domains
    densitytensor = DensityTensor(data, domains, measurement_names)

    return densitytensor
end

"""
    names(n::NamedArray, dimname::Union{String,Symbol})

Extend the names function from NamedArray to get the names given the axis name rather than
the axis number.
"""
function Base.names(n::NamedArray, dimname::Union{String,Symbol})
    index = findfirst(dimnames(n) .== dimname)
    if isnothing(index)
        KeyError("Dimention $dimname not found in array.")
    end
    return names(n, index)
end

# Getters for useful quantities
getmeasurements(D::DensityTensor) = names(D, "measurement")

function _getmeasurementindex(D::DensityTensor, measurement::String)
    index = findfirst(names(D, "measurement") .== measurement)
    if isnothing(index)
        KeyError("Measurment $measurement not found in DensityTensor.")
    end
    return index
end

"""
    getdomain(D::DensityTensor, measurement::String)
    getdomain(D::DensityTensor, j::Integer)

Gets the domain for the `measurement` density, the locations where the density was
sampled.
See [`getdomains`](@ref).
"""
getdomain(D::DensityTensor, measurement::String) = getdomain(D, _getmeasurementindex(D, measurement))
getdomain(D::DensityTensor, j::Integer) = getdomains(D)[j]

"""
    getdomains(D::DensityTensor{T})::Vector{Vector{T}}

Gets the domain for every measurement's density, the locations where each density was
sampled.
See [`getdomains`](@ref).
"""
function getdomains(D::DensityTensor)
    return collect.(zip(names(D, "density")...))
end

"""
    getsource(D::DensityTensor, i::Integer)
    getsink(D::DensityTensor, i::Integer)

Gets source/sink i from D. See [`eachsource`](@ref).
"""
getsource(D::DensityTensor, i::Integer) = D[i, :, :] # TODO see if @view is better

"""Alias for [`getsource`](@ref)"""
const getsink = getsource

"""
    getstepsizes(D::DensityTensor)

Gets the step sizes used for each domain. See [`getdomains`](@ref).
"""
getstepsizes(D::DensityTensor) = [x[begin+1] - x[begin] for x in getdomains(D)]

"""
    getsourcename(D::DensityTensor)

Gets the name for the grouping of measurements. Usually `"Sink"` or `"Source"`.
See [`getsourcename`](@ref).
"""
getsourcename(D::DensityTensor) = dimnames(D)[1]

"""
    getsourcenames(D::DensityTensor)

Gets the list of all sources' names. For example, `["Sink 1", "Sink 2", ...]`.
See [`getsourcenames`](@ref).
"""
getsourcenames(D::DensityTensor) = names(D, 1)#["$(getsourcename(D)) $s" for s in ]

# Setters
"""
    setsourcename!(D::DensityTensor, name::String)

Sets the name of the source used by [`getsourcename`](@ref) and [`getsourcenames`](@ref).
"""
setsourcename!(D::DensityTensor, name::String) = setdimnames!(D::NamedArray, name, 1)

# Other manipulators
"""
    normalize_density_sums!(D::DensityTensor)

Rescales the densities so that the sum of the density samples is 1.

This is in constrast to the usualy normalization for density functions where the area of the
density curve is 1. In the case of an evenly sampled density, this area is
`sum(density_samples)*step_size`.

Use [`normalize_density_sums`](@ref) to avoid mutation.
"""
function normalize_density_sums!(D::DensityTensor)
    for (slice, x) in zip(eachmeasurement(D), getstepsizes(D))
        slice .*= x
    end
end

"""See [`normalize_density_sums!`](@ref)"""
function normalize_density_sums(D::DensityTensor)
    D_copy = deepcopy(D)
    normalize_density_sums!(D_copy)
    return D_copy
end

# Iterators
"""
    eachdensity(D::DensityTensor)
    eachdensity(D::DensityTensor, measurement::String)

Iterates D over each density vector. These are the 3 fibers of D.
If a measurement is given, iterates over the densities for that measurement.
"""
eachdensity(D::DensityTensor) = eachslice(D, dims=(1,2))
eachdensity(D::DensityTensor, measurement::String) = eachrow(D[:, measurement, :])

"""
    eachmeasurement(D::DensityTensor)

Iterates D over each measurement slice. These are the lateral slices.
"""
eachmeasurement(D::DensityTensor) = eachslice(D, dims=2)

"""
    eachsource(D::DensityTensor)
    eachsink(D::DensityTensor)

Iterates D over each source/sink slice. These are the horizontal slices.
See [`getsource`](@ref).
"""
eachsource(D::DensityTensor) = eachslice(D, dims=1)

"""Alias for [`eachsource`](@ref)."""
const eachsink = eachsource
