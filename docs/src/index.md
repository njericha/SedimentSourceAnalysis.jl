# Sediment Source Analysis
<!--
```@meta
DocTestSetup = quote
    using SedimentAnalysis
end
```
-->

```@contents
Depth = 3
```
## Importing
For now, clone the repo and install with `julia> ] add SedimentAnalysis`.
To import, type `using SedimentAnalysis` at the top of your file to load both submodules.
If only one of the modules is desired, type `using SedimentAnalysis.XXX`.
The modules are built to be independent of each other so that (eventually) the MTF could be moved to an separate package altogether.

## Submodules
The two main submodules are MTF (**M**atrix **T**ensor **F**actorization) and SedimentTools.

### MTF
Defines the main factorization function [`nnmtf`](@ref) and related mathematical functions. See the full documentation here [Matrix Tensor Factorization](@ref).

### SedimentTools
Holds various types at the [`Grain`], and [`Sink`] level, importing ([`read_raw_data`]) and processing data ([`make_densities`]) functions, and additional methods of some [Plots.jl](https://docs.juliaplots.org/stable/) functions for visualization with these custom types.

## Index

```@index
```
