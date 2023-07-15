#= Single module that pulls from the
factorization and data handling modules =#

module SedimentAnalysis #TODO check if these files should be included and reexported instead

include("MTF/MTF.jl")
include("SedimentTools/SedimentTools.jl")
using .MTF
using .SedimentTools

export * # reexport everything

end # SedimentAnalysis
