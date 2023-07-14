#=
main script for factorizing the data at the grain level
=#
using Plots
using StatsPlots
using Statistics
using NMF

# Import data and list of elements
include("./dataimport.jl")
raw_data_path(s::String) = join(["rawdata/Lee et al 2021 " s ".xlsx"],"")
elements = ["Y", "Nb", "Ce", "Pr", "Nd", "Sm", "Eu_Euprime", "U", "age"]
grain_data, rock_index = read_xlsx_data(raw_data_path, elements);
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
rank = 5
results = nnmf(G, rank; alg=:cd, maxiter=400, tol=1e-3)

# Extract coefficient and factor matricies
C, F = results.W, results.H

#= Rescale the matricies to better interperate data =#
# Normalize factors in latent space
factor_norms = Diagonal(norm.(eachrow(F)))
C = C * factor_norms
F = factor_norms^(-1) * F
# Ensure row sum of coefficients is 1
substance_amount = Diagonal(dropdims(sum(C,dims=2),dims=2))
coefficients = substance_amount^(-1) * C
# Rescale factors from latent space to ppm
factors = F * S

# Evaluate performance
objective = 0.5 * norm(G - C*F)^2
rel_error_obj = norm(G - C*F)/norm(G)
rel_error_data = norm(grain_data - substance_amount*coefficients*factors)/norm(grain_data)

# Visualizing factors
heatmap(F,
    yticks=(1:rank, ["factor $i" for i ∈ 1:rank]),
    xticks=(1:length(elements), elements),
    title = "Latent Factors in Grain Data")

heatmap(coefficients,
    yticks=(start_index, ["rock $i" for i ∈ 1:length(rock_index)]),
    xticks=(1:rank,["factor $i" for i ∈ 1:rank]),
    title="Proportion of Factors in Grain Samples",
    size=(500,600))

import Base: split
split(v::Vector{T}, starts_and_ends) where T = [v[i:j] for (i,j) ∈ starts_and_ends]

# Split coefficient data by rock for each factor
for j ∈ 1:rank
    fs = split(coefficients[:,j], rock_index)
    p = violin(fs,
        xticks=1:length(rock_index),
        ylabel="Coefficient",
        xlabel="Rock #",
        title="Proportion of Factor $j in Each Rock",
        legend= false)
    display(p)
end