using Documenter

push!(LOAD_PATH,"../src/")
#using Pkg
using BlockTensorDecomposition
using SedimentSourceAnalysis
#using SedimentSourceAnalysis.MTF
using SedimentSourceAnalysis.SedimentTools
using NamedArrays
using Plots

DocMeta.setdocmeta!(
    SedimentSourceAnalysis,
    :DocTestSetup,
    :(using SedimentSourceAnalysis;
    # using SedimentSourceAnalysis.SedimentTools;
    using Pkg; Pkg.add(url="https://github.com/MPF-Optimization-Laboratory/BlockTensorDecomposition.jl.git");
    using BlockTensorDecomposition: make_densities; using NamedArrays: NamedArray; using Plots: heatmap);
    recursive=true
)

makedocs(
    sitename="Sediment Source Analysis",
    modules = [SedimentSourceAnalysis, SedimentTools], #MTF
    checkdocs = :exports,
    repo = Documenter.Remotes.GitHub("njericha", "SedimentSourceAnalysis.jl")
)

deploydocs( # TODO either use SedimentSourceAnalysis or SedimentSourceAnalysis
    repo = "github.com/njericha/SedimentSourceAnalysis.jl.git",
)
