struct DistributionData{T}
    name::String
    scale::AbstractVector{T}
    data::AbstractMatrix{T}
end

struct Distributions{T} <: AbstractVector{DistributionData{T}}
    names::Vector{String}
    scales::Vector{AbstractVector{T}}
    data::Vector{AbstractMatrix{T}}
end

function Distributions(v::AbstractVector{DistributionData})
    names = [d.name for d ∈ v]
    scales = [d.scale for d ∈ v]
    data = [d.data for d ∈ v]
    return Distributions(names, scales, data)
end

import Base: size, getindex
Base.size(D::Distributions) = (length(D.names),)
Base.getindex(D::Distributions, i::Int) = DistributionData(D.names[i], D.scales[i], D.data[i])

"""
Converts the Distributions to a plain Array{T,3} type
"""
to_tensor(D::Distributions) = cat(D.data..., dims=3)


############
#=
struct Grain{T} <: AbstractVector{T} where T <: Number
    ordereddict::OrderedDict{String, T}
end

# Interface
# allow for alternate backend like NamedTuple, NamedArray, etc.
import OrderedCollections: OrderedDict
OrderedCollections.OrderedDict{String,T}(g::Grain{T}) where T <: Number = g.ordereddict

import Base: size, getindex, keys, values, collect, size, getindex, iterate, convert, setindex!
Base.Vector{T}(g::Grain{T}) where T <: Number = values(g)

"""List of names for each measurement"""
Base.keys(g::Grain) = keys(OrderedDict(g))

"""A plain Vector{T} with the values for each measurment"""
Base.values(g::Grain) = values(OrderedDict(g))

Base.collect(g::Grain) = collect(OrderedDict(g))
Base.size(g::Grain) =(length(Vector(g)),)

"""gets the i'th (name, value) pair in the grain"""
Base.getindex(g::Grain, key::String) = OrderedDict(g)[key]
Base.setindex!(g::Grain{T}, value::T, key::String) = setindex!(OrderedDict(g), value, key)
=#

#########
# TODO use NamedArray backend
# Grains
Grain{T <: Number} = OrderedDict{String, T}
import Base: Vector

"""turns a grain into a regular Vector type"""
Base.Vector{T}(g::Grain{T}) where T <: Number = collect(values(g))
Vector(g::Grain{T}) where T <: Number = Vector{T}(g)
getnames(g::Grain) = collect(keys(g))
getvalues(g::Grain) = collect(values(g))
Grain(iter) = OrderedDict(iter)

# Sinks
# Want to ensure the grains have the same names and are in the same order
Sink{T <: Number} = Vector{Grain{T}} # Using vector and not a set to preserve order
#Rock = Sink # TODO Make Alias and check it works

getnames(s::Sink) = iszero(length(s)) ? String[] : get_names(s[1])
getindex(s::Sink, key::String) = [g[key] for g ∈ s]

"""
    Sink(grain1, grain2, ...)
    Sink([grain1, grain2, ...])

Collects a list of Grains into a Rock/Sink.

Ensures all Grains have the same names and are in the same order.
"""
function Sink(vec_of_grains::AbstractVector{Grain}) #want each element to be a grain
    @assert allequal(getnames.(vec_of_grains))
    return collect(vec_of_grains)::Sink
end
Sink(vec_of_grains::AbstractVector{Grain}...) = Sink(vec_of_grains)

#=
Alternatively use:
struct Sink{T} <: AbstractVector{Grain{T}} where T <: Number
    names::AbstractVector{String}
    grains::AbstractVector{Grain{T}}
end
get_names(s::Sink) = s.names
function Sink(vec_of_grains::AbstractVector{Grain}) #want each element to be a grain
    @assert allequal(get_names.(vec_of_grains))
    return Sink(get_names(vec_of_grains[1]), vec_of_grains)
end
=#
# TODO Should this really be part of the Sink struct where the distributions
# are auto calculated? May not want to because of it would hide the computation.
struct DistributionSink{T} <: AbstractArray{T,3}
    names::AbstractVector{String}
    distributions::AbstractArray{T,3}
    bandwidths::AbstractVector{<:Real}
    scales::AbstractVector{AbstractVector{T}}
end

getnames(d::DistributionSink) = d.names
getdistributions(d::DistributionSink) = d.distributions
getbandwidths(d::DistributionSink) = d.distributions
getscales(d::DistributionSink) = d.scales
getindex(d::DistributionSink, args...) = getindex(getdistributions(d), args...)
size(d::DistributionSink) = size(getdistributions(d))
to_tensor(d::DistributionSink) = collect(getdistributions(d))

"""
    DistributionSink(s::Sink, steps::Integer; kwargs...)

Construct a DistributionSink from a sink.

See [`make_distributions`](@ref) for allowed keyword arguments.
"""
function DistributionSink(
    s::Sink,
    steps::Integer;
    kwargs...
)
    # Check input is in the correct range
    (0 < inner_percentile <= 100) ||
        ArgumentError("inner_percentile must be between 0 and 100, got $inner_percentile")

    names = getnames(s)
    distributions, bandwidths, scales = make_distributions(s, steps; kwargs...)
    return DistributionSink(names, distributions, bandwidths, scales)
end

"""Shortcut to get distributions streight from a Sink"""
function getdistributions(s::Sink, steps; kwargs...)
    return getdistributions(DistributionSink(s, steps; kwargs...))
end
