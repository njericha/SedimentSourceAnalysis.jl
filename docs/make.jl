using Documenter

#push!(LOAD_PATH,"../src/")
using SedimentAnalysis
#using SedimentAnalysis.MTF
using SedimentAnalysis.SedimentTools
using MatrixTensorFactor
using NamedArrays
using Plots

DocMeta.setdocmeta!(
    SedimentAnalysis,
    :DocTestSetup,
    :(using SedimentAnalysis; using MatrixTensorFactor; using NamedArrays: NamedArray; using Plots: heatmap);
    recursive=true
)

makedocs(
    sitename="Sediment Source Analysis",
    modules = [SedimentAnalysis, SedimentTools], #MTF
)

deploydocs( # TODO either use SedimentSourceAnalysis or SedimentAnalysis
    repo = "github.com/njericha/Sediment-Source-Analysis.jl.git",
)
