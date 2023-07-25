#= Single module that combines the
factorization and data handling modules =#

module SedimentAnalysis #TODO check if these files should be included or imported instead

include("MTF/MTF.jl")
include("SedimentTools/SedimentTools.jl")

using .MTF
using .SedimentTools

export combined_norm, nnmtf, plot_factors, rel_error # matrixtensorfactorization

export Grain, DensityTensor, Rock, Sink, Source # Types
export array, getdomain, getdomains, getsourcename, getsourcenames, getmeasurements
export nammedarray, sink, source # Getters
export normalize_density_sums!, normalize_density_sums, setsourcename! # Setters
export eachdensity, eachmeasurement, eachsink, eachsource # Iterators

export read_raw_data

export DEFAULT_ALPHA, DEFAULT_N_SAMPLES # Constants
export default_bandwidth, make_densities, standardize_KDEs # Functions

export estimate_which_source, match_sources!

export measurement_heatmaps, plot_densities, source_heatmaps

end # module SedimentAnalysis
