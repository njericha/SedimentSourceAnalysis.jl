# Sediment Source Analysis

```@contents
Depth = 3
```
## How to use

### In Browser
1. Go to https://github.dev/njericha/Sediment-Source-Analysis.jl

Note you can start a terminal with `Ctrl+Shift+C` (Windows) or `Cmd+Shift+C` (Mac).

**OR**
### On your own device
1. Clone the repo at https://github.com/njericha/Sediment-Source-Analysis.jl 

### Importing the package
2. Install with `julia> ] add SedimentAnalysis`.
3. To import into a file or in the REPL, type `using SedimentAnalysis` load both submodules (`MTF` and `SedimentTools`).
4. If only one of the modules is desired, type `using SedimentAnalysis.XXX`.

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
