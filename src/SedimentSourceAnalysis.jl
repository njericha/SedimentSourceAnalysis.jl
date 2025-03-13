#= Single module that combines the
factorization and data handling modules =#

module SedimentSourceAnalysis #TODO check if these files should be included or imported instead

#include("MTF/MTF.jl")

include("SedimentTools/SedimentTools.jl")

#using .MTF
using .SedimentTools

# re-exported from BlockTensorDecomposition
export combined_norm, dist_to_Ncone, nnmtf, plot_factors, rel_error, mean_rel_error
export d_dx, d2_dx2, curvature, standard_curvature # Approximations
#export make_densities, make_densities2d # method extention

export Grain, DensityTensor, Rock, Sink, Source # Types
export array, getdomain, getdomains, getsourcename, getsourcenames, getmeasurements
export getstepsizes, namedarray, getsink, getsource # Getters
export normalize_density_sums!, normalize_density_sums, setsourcename! # Setters
export eachdensity, eachmeasurement, eachsink, eachsource # Iterators

export read_raw_data

#export DEFAULT_ALPHA, DEFAULT_N_SAMPLES # Constants
#export default_bandwidth, make_densities, standardize_KDEs # Functions

export estimate_which_source, estimate_which_2d_source, estimate_which_nd_source, label_accuracy, match_sources!, confidence_score

export measurement_heatmaps, plot_densities, source_heatmaps, plot_convergence, plot_source_index

end # module SedimentSourceAnalysis
