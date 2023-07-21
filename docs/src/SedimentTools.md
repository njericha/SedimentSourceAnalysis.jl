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
measurments
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
eachmeasurment
eachsink
eachsource
```

### Density Estimation
```@docs
default_bandwidth
make_densities
standardize_KDEs
```

### Visualizers
```@docs
measurment_heatmaps
source_heatmaps
```

## Index

```@index
```
