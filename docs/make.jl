using Documenter

#push!(LOAD_PATH,"../src/")
using SedimentAnalysis
using SedimentAnalysis.MTF
using SedimentAnalysis.SedimentTools
using NamedArrays

DocMeta.setdocmeta!(
    SedimentAnalysis,
    :DocTestSetup,
    :(using SedimentAnalysis; using NamedArrays: NamedArray);
    recursive=true
)

makedocs(
    sitename="Sediment Source Analysis",
    modules = [SedimentAnalysis, MTF, SedimentTools],
)

deploydocs( # TODO either use SedimentSourceAnalysis or SedimentAnalysis
    repo = "github.com/njericha/Sediment-Source-Analysis.jl.git",
)
