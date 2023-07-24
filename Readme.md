# Sediment Source Analysis

<!---
comment text
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://njericha.github.io/SedimentAnalysis.jl/stable)



GitHub Actions : [![Build Status](https://github.com/JuliaLang/Example.jl/workflows/CI/badge.svg)](https://github.com/njericha/Sediment-Source-Analysis/actions?query=workflow%3ACI+branch%3Amaster)

AppVeyor: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaLang/Example.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tkelman/example-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/JuliaLang/Example.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaLang/Example.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLang/Example.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLang/Example.jl?branch=master)
--->

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://njericha.github.io/Sediment-Source-Analysis.jl/dev/)

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