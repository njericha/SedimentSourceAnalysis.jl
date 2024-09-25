# Sediment Source Analysis

<!---
comment text
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://njericha.github.io/SedimentSourceAnalysis.jl/stable)



GitHub Actions : [![Build Status](https://github.com/JuliaLang/Example.jl/workflows/CI/badge.svg)](https://github.com/njericha/SedimentSourceAnalysis/actions?query=workflow%3ACI+branch%3Amaster)

AppVeyor: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaLang/Example.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tkelman/example-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/JuliaLang/Example.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaLang/Example.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLang/Example.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLang/Example.jl?branch=master)
--->

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://njericha.github.io/SedimentSourceAnalysis.jl/dev/)

## Summary of Features
- Import grain data from Excel files
- Perform kernel density estimation on grain data by sinks and features
- Decompose sink distributions into weighted combination of sources
- Label grains based on learned sources
- Plot kernel densities, sink, and source distributions

## Repo organization

Our main numerical experiment can be found in `knownsources median bandwidth.jl` under the `examples` folder.

The data used in our main numerical experiment `knownsources median bandwidth.jl` can be found [here](https://github.com/njericha/SedimentSourceAnalysis.jl/tree/main/data/sundell2022) under `data/sundell2022`. The exact grains selected and number of grains from each sink can be found in `data/20sinks from 3Sources from Sundell et al 2022.xlsx`.

The data handling and backend functions can be found in `src/SedimentTools`.

Our main decomposition algorithm can be found in the more general repo [MatrixTensorFactor.jl](https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl).


## How to run the code

### In Browser
1. Go to https://github.com/njericha/SedimentSourceAnalysis.jl
2. Click "<> Code" and press "+" to "Create a codespace on main". It make take a few moments to set up.
3. Open the command palett with `Ctrl+Shift+P` (Windows) or `Cmd+Shift+P` (Mac)
4. Enter `>Julia: Start REPL`
5. In the REPL, resolve any dependency issues with `pkg> resolve` and `pkg> instantiate` (use `julia> ]` to get to the package manager). It may take a few minutes to download dependencies.

Run one of the example files by opening the file and pressing the triangular "run" button, or `>Julia: Execute active File in REPL`.

**OR**
### On your own device
1. Clone the repo at https://github.com/njericha/SedimentSourceAnalysis.jl
2. Navigate to the root of the repository in a terminal and run `julia`
3. Activate the project with `pkg> activate .` (use `julia> ]` to get to the package manager)
4. resolve any dependency issues with `pkg> resolve`

### Importing the package
Type `julia> using SedimentSourceAnalysis` load the package, or `using SedimentSourceAnalysis.SedimentTools` to load the submodule directly.

## Examples
See the `examples` folder for the following files.

`knownsources median bandwidth.jl`: Uses data from Sundel et al where we know the sources of each Grain. Use this to see how well the factorization performs with realistic data.

`knownsources.jl`: Similar to `knownsources median bandwidth.jl`, but uses the bandwidth from the first sink.

`measurementcorrelation.jl`: Checks the validity of representing the grain distributions as a product distribution. In particular, we would like the measurements to be independent.

`unknownsources.jl`: Uses data from Lee et al where we don't have a ground truth. Showcases how the method would be used in practice.

`randomtensor.jl`: Factorizes a random 50x50x50 tensor. See how the factorization performs in theory when a perfect factorization exists.

## Submodules
The main submodule of this repo is SedimentTools. The submodule MTF has been moved to a separate repo [MatrixTensorFactor](https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl).

### SedimentTools
Holds various types at the [`Grain`], and [`Sink`] level, importing ([`read_raw_data`]) and processing data ([`make_densities`]) functions, and additional methods of some [Plots.jl](https://docs.juliaplots.org/stable/) functions for visualization with these custom types.

### MatrixTensorFactor
Defines the main factorization function [`nnmtf`](@ref) and related mathematical functions. See the repo here [MatrixTensorFactor.jl](https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl).

# Citation

If you find this repo helpful, please cite the associated paper:

```
@misc{graham_tracing_2024,
  title = {Tracing {Sedimentary} {Origins} in {Multivariate} {Geochronology} via {Constrained} {Tensor} {Factorization}},
  url = {https://friedlander.io/publications/2024-sediment-source-analysis/},
  author = {Graham, Naomi and Richardson, Nicholas and Friedlander, Michael P. and Saylor, Joel},
  month = may,
  year = {2024},
}
```

or feel free to reach out to us with an email to njericha at math.ubc.ca.
