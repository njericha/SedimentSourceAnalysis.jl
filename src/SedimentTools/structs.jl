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
