"""
Import, store, and visualize sediment data.
"""
module SedimentTools

using XLSX # TODO annotate exactly what is being used from each imported package

#export ...

include("importing.jl")
include("structs.jl")
include("densityestimation.jl")
include("viz.jl")
end
