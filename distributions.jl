#=
distributions.jl

Defines the custom struct to keep track of the distribution data and soem interfacing
=#

struct DistributionData{T}
    name::String
    scale::AbstractVector{T}
    data::AbstractMatrix{T}
end

struct Distributions{T} <: AbstractVector{DistributionData{T}}
    names::Vector{String}
    scales::Vector{V} where V <: AbstractVector{T}
    data::Vector{M} where M <: AbstractMatrix{T}
end

function Distributions(v::AbstractVector{DistributionData})
    names = [d.name for d ∈ v]
    scales = [d.scale for d ∈ v]
    data = [d.data for d ∈ v]
    return Distributions(names, scales, data)
end

import Base: size, getindex, convert
Base.size(D::Distributions) = (length(D.names),)
Base.getindex(D::Distributions, i::Int) = DistributionData(D.names[i], D.scales[i], D.data[i])

# Add method to `convert` so I can read strings as Floats
# This should not be used on arbitrary strings and only on numbers surrounded
# by double quotes. This is nessesary because some excel entries have mistakenly
# read the number as a string.
Base.convert(::Type{T}, s::String) where T <: Number = parse(T, s)

"""
Converts the Distributions to a regular Array{T,3} type
"""
to_tensor(D::Distributions) = cat(D.data..., dims=3)

"""
Strip units like "_ppm" or "_Ma" from the long name titles
"""
function short_name(s::String)
    if s == "Dy_Yb"
        return "Dy-Yb" # leave the Dy Yb element combination together
    else
        return split(s,"_")[1]
    end    
end