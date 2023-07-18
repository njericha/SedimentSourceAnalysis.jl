#= Single module that combines the
factorization and data handling modules =#

module SedimentAnalysis #TODO check if these files should be included and reexported instead
using Reexport

include("MTF/MTF.jl")
include("SedimentTools/SedimentTools.jl")
#@reexport using .MTF
#@reexport using .SedimentTools
using .MTF
export nnmft, combined_norm

using .SedimentTools

end # SedimentAnalysis
