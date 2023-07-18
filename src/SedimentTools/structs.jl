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

struct Grain{T} <: AbstractVector{T}
    ordereddict::OrderedDict{String, T}
end

# Interface
# allow for alternate backend like NamedTuple, NamedArray, etc.
OrderedDict{String,T}(g::Grain{T}) where T <: Number = g.ordereddict

import Base: size, getindex, keys, values, collect, size, getindex, iterate, convert
Base.Vector{T}(g::Grain{T}) where T <: Number = values(g)

"""List of names for each measurement"""
Base.keys(g::Grain) = keys(OrderedDict(g))

"""A plain Vector{T} with the values for each measurment"""
Base.values(g::Grain) = values(OrderedDict(g))

Base.collect(g::Grain) = collect(OrderedDict(g))
Base.size(g::Grain) =(length(Vector(g)),)

"""gets the i'th (name, value) pair in the grain"""
Base.getindex(g::Grain, key::String) = OrderedDict(g)[key]
