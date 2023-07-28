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

## How setup the environment

### In Browser
1. Go to TODO UPDATE: https://github.com/codespaces/njericha-expert-space-goggles-96w695j9q5g274r7
2. Open the command palett with `Ctrl+Shift+P` (Windows) or `Cmd+Shift+P` (Mac)
3. Enter `>Julia: Start REPL`
4. In the REPL, resolve any dependency issues with `pkg> resolve` (use `julia> ]` to get to the package manager)

Run one of the example files by opening the file and pressing the triangular "run" button, or `>Julia: Execute active File in REPL`.

**OR**
### On your own device
1. Clone the repo at https://github.com/njericha/Sediment-Source-Analysis.jl
2. Navigate to the root of the repository in a terminal and run `julia`
3. Activate the project with `pkg> activate .` (use `julia> ]` to get to the package manager)
4. resolve any dependency issues with `pkg> resolve`

### Importing the package
Type `julia> using SedimentAnalysis` load both submodules (`MTF` and `SedimentTools`), or if only one of the modules is desired, type `using SedimentAnalysis.XXX`.

The modules are built to be independent of each other so that (eventually) the MTF could be moved to an separate package altogether.

## Examples
`knownsources.jl`: Uses data from Sundel et al where we know the sources of each Grain. Use this to see how well the factorization performs with realistic data.
`unknownsources.jl`: Uses data from Lee et al where we don't have a ground truth. Showcases how the method would be used in practice.
`randomtensor`: Factorizes a random 50x50x50 tensor. See how the factorization performs in theory when a perfect factorization exists.

## Submodules
The two main submodules are MTF (**M**atrix **T**ensor **F**actorization) and SedimentTools.

### MTF
Defines the main factorization function [`nnmtf`](@ref) and related mathematical functions. See the full documentation here [Matrix Tensor Factorization](@ref).

### SedimentTools
Holds various types at the [`Grain`], and [`Sink`] level, importing ([`read_raw_data`]) and processing data ([`make_densities`]) functions, and additional methods of some [Plots.jl](https://docs.juliaplots.org/stable/) functions for visualization with these custom types.
