using Documenter

#push!(LOAD_PATH,"../src/")
using Pkg
Pkg.add(url="https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl.git")
using MatrixTensorFactor
using SedimentAnalysis
#using SedimentAnalysis.MTF
using SedimentAnalysis.SedimentTools
using NamedArrays
using Plots

DocMeta.setdocmeta!(
    SedimentAnalysis,
    :DocTestSetup,
    :(using SedimentAnalysis;
    using Pkg; Pkg.add(url="https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl.git");
    using MatrixTensorFactor: make_densities; using NamedArrays: NamedArray; using Plots: heatmap);
    recursive=true
)

makedocs(
    sitename="Sediment Source Analysis",
    modules = [SedimentAnalysis, SedimentTools], #MTF
    checkdocs = :exports,
)

deploydocs( # TODO either use SedimentSourceAnalysis or SedimentAnalysis
    repo = "github.com/njericha/Sediment-Source-Analysis.jl.git",
)
