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

# TODO use NamedArray backend

#########

# Grains
#=
Grain{T <: Number} = OrderedDict{String, T}
import Base: Vector

"""turns a grain into a regular Vector type"""
Base.Vector{T}(g::Grain{T}) where T <: Number = collect(values(g))
Vector(g::Grain{T}) where T <: Number = Vector{T}(g)
names(g::Grain) = collect(keys(g))
values(g::Grain) = collect(values(g::OrderedDict))
Grain(iter) = OrderedDict(iter)

# Sinks
# Want to ensure the grains have the same names and are in the same order
Sink{T <: Number} = Vector{Grain{T}} # Using vector and not a set to preserve order
#Rock = Sink # TODO Make Alias and check it works

names(s::Sink) = iszero(length(s)) ? String[] : names(s[1])
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
    names::Vector{String}
    distributions::Vector{UnivariateKDE}
    bandwidths::Vector{<:Real}
    scales::Vector{AbstractVector{T}}
    sampled_distributions::AbstractArray{T,3}
end

names(d::DistributionSink) = d.names
distributions(d::DistributionSink) = d.sampled_distributions
bandwidths(d::DistributionSink) = d.bandwidths
scales(d::DistributionSink) = d.scales
getindex(d::DistributionSink, args...) = getindex(distributions(d), args...)
size(d::DistributionSink) = size(distributions(d))

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
    sampled_distributions, distributions, bandwidths, scales = make_distributions(s, steps; kwargs...)
    return DistributionSink(names(s), distributions, bandwidths, scales, sampled_distributions)
end

"""Shortcut to get distributions streight from a Sink"""
function distributions(s::Sink, steps; kwargs...)
    return distributions(DistributionSink(s, steps; kwargs...))
end

to_tensor(sinks::AbstractVector{Sink}) = cat(collect.(distributions(sinks))..., dims=1)
=#

#########
# Grain #
#########

Grain{T <: Number} = NamedVector{T}
measurments(g::Grain) = names(g, 1) # names(g) from NamedArray returns Vector{Vector{T}}

#################
# Sinks / Rocks #
#################

Sink{T <: Number} = Vector{Grain{T}} # Using vector and not a set to preserve order
Rock = Sink

measurments(s::Sink) = iszero(length(s)) ? String[] : measurments(s[1])
getindex(s::Sink, key::String) = (g[key] for g ∈ s)

"""Iterator for a list of values of each measurement"""
eachmeasurment(s::Sink) = (s[m] for m in measurments(s))

"""
    Sink(grain1, grain2, ...)
    Sink([grain1, grain2, ...])

Collects a list of Grains into a Rock/Sink.

Ensures all Grains have the same names and are in the same order.
"""
function Sink(vec_of_grains::AbstractVector{Grain}) #want each element to be a grain
    @assert allequal(measurments.(vec_of_grains))
    return collect(vec_of_grains)::Sink
end
Sink(vec_of_grains::AbstractVector{Grain}...) = Sink(vec_of_grains)

####################
# DistributionSink #
####################

DistributionSinks{T <: Number} = NamedArray{T, 3}

measurments(d::DistributionSinks) = names(d, "measurments")
getsink(d::DistributionSinks, i::Integer) = d[i, :, :] # TODO see if @view is better
getsource = getsink

"""Methods to get the names given the axis name rather than the axis number"""
function Base.names(n::NamedArray, dimname::Union{String,Symbol})
    return names(n, findfirst(dimnames(n) .== dimname))
end
Base.names(n::NamedArray, dimname::Name) = names(n, findfirst(dimnames(n) .== dimname.names))

to_tensor(sinks::DistributionSinks) = sinks.array

function DistributionSinks(sinks::AbstractVector{Sinks}, KDEs::AbstractVector{AbstractVector{UnivariateKDE}})
    # Argument Handeling
    allequal(measurments.(sinks)) ||
        ArgumentError("All sinks must have the same measurements in the same order.")
    length(sinks) == length(KDEs) ||
        ArgumentError("Must be the same number of sinks as there are lists of KDEs.")
    measurment_names = measurments(sinks[begin])
    length(measurment_names) == length(KDEs[begin]) ||
        ArgumentError("Must be the same number of measurements as there are KDEs for each sink.")

        # TODO make this line more legible, possible by naming the KDEs to begin with
    # Magic line to take the KDEs into an order-3 tensor
    data = permutedims(cat(cat(map.(k -> k.density, KDEs)..., dims=2)..., dims=3), [3,2,1])

    # Confirm all the dimentions are in the right order
    n_density_samples = length((KDEs[begin][begin]).x)
    @assert size(data) == (length(sinks), length(measurment_names), n_density_samples)

    # Wrap in a NamedArray
    dsinks = NamedArray(data, dimnames=("sink", "measurment", "density"))::DistributionSinks
    setnames!(dsinks, measurment_names, 2)

    return dsinks
end
