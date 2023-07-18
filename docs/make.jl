using Documenter

#push!(LOAD_PATH,"../src/")
using SedimentAnalysis
using SedimentAnalysis.MTF
using SedimentAnalysis.SedimentTools

# DocMeta.setdocmeta!(
#     SedimentAnalysis,
#     :DocTestSetup,
#     :(using .MTF);
#     recursive=true
# )

makedocs(
    sitename="Sediment Source Analysis",
    #modules = [SedimentAnalysis], # could use [MTF, SedimentTools] ?
    modules = [SedimentAnalysis, MTF, SedimentTools], # could use  ?
)

deploydocs(
    repo = "github.com/njericha/Sediment-Source-Analysis.jl.git",
)
