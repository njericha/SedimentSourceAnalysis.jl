#=
main script for factorizing the data at the grain level
    for the  known sources
=#
using Plots
using StatsPlots
using Statistics
using NMF

# Import data and list of elements
include("./dataimport.jl")
elements = ["ages", "Eu_anomaly", "Ti_temp", "LREE_HREE", "Dy_Yb", "(Ce_Nd)_Y"]#"Th_U"
data_path = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
grain_data, rock_index = import_grain_data(data_path, elements);
start_index = [idx for (idx, _) ∈ rock_index]
end_index = [idx for (_, idx) ∈ rock_index]

# Normalize each column so they are weighted equaly
# This ensures the factorization algorithm doesn't
# favour fitting one element over another
column_norm(A) = norm.(eachcol(A))
column_median(A) = dropdims(median(A, dims = 1),dims=1)
column_max(A) = maximum.(eachcol(A))
column_scale(A; δ=1e-10) = Diagonal(column_median(A) .+ δ) # δ ensures safe inversion

S = column_scale(grain_data)
G = (grain_data * S^(-1))
G = convert(Matrix{Float64},G) # ensure type stability from lingering Matrix{Any} types

# Add small number for stability
# Some methods cannot have values identicaly zero
G .+= 1.0e-10

# Perform the low-rank decomposition
rank = 3
include("./coorddecent.jl")
C, F, error, norm_grad = binarycoorddecent(G, rank, maxiter=300)
#factors = F * S

# Evaluate performance
objective = 0.5 * norm(G - C*F)^2
rel_error_obj = norm(G - C*F)/norm(G)
#rel_error_data = norm(grain_data - substance_amount*coefficients*factors)/norm(grain_data)

# Visualizing factors
heatmap(F,
    yticks=(1:rank, ["source $i" for i ∈ 1:rank]),
    xticks=(1:length(elements), elements),
    title = "Latent Factors in Grain Data")

heatmap(C,
    yticks=(start_index, ["sink $i" for i ∈ 1:length(rock_index)]),
    xticks=(1:rank,["source $i" for i ∈ 1:rank]),
    title="Proportion of Factors in Grain Samples",
    size=(500,600))

import Base: split
split(v::Vector{T}, starts_and_ends) where T = [v[i:j] for (i,j) ∈ starts_and_ends]

# Split coefficient data by rock for each factor
for j ∈ 1:rank
    fs = split(C[:,j], rock_index)
    p = violin(fs,
        xticks=1:length(rock_index),
        ylabel="Coefficient",
        xlabel="Rock #",
        title="Proportion of Factor $j in Each Rock",
        legend= false)
    display(p)
end

# Calculuate C_true
#
