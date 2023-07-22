# Sediment Analysis Tools

```@contents
depth = 3
```

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
heatmap(::NamedArray)
```

## Functions

### Importers
```@docs
read_raw_data
```

### Getters
```@docs
array
domain
domains
getsourcename
measurements
nammedarray
sink
source
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
```

## Index

```@index
```
