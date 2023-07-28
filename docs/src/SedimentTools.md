# Sediment Analysis Tools

## Types
```@docs
Grain
DensityTensor
Rock
Sink
Source
```

## Constants
```@docs
DEFAULT_ALPHA
DEFAULT_N_SAMPLES
```

## Method Extentions
```@docs
getindex(::Sink, ::String)
names(::NamedArray, ::Union{String,Symbol})
heatmap(::NamedMatrix)
```

## Functions

### Importers
```@docs
read_raw_data
```

### Getters
```@docs
array
getdomain
getdomains
getmeasurements
getsourcename
getsourcenames
getstepsizes
namedarray
getsink
getsource
```

### Setters and Manipulators
```@docs
normalize_density_sums!
normalize_density_sums
setsourcename!
```

### Iterators
```@docs
eachdensity
eachmeasurement
eachsink
eachsource
```

### Density and Source Estimation
```@docs
estimate_which_source
default_bandwidth
make_densities
match_sources!
standardize_KDEs
```

### Visualizers
```@docs
measurement_heatmaps
plot_densities
source_heatmaps
plot_convergence
plot_source_index
```

## Index

```@index
```
