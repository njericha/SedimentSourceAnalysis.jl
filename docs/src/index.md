# Sediment Source Analysis

## Summary of Features
- Import grain data from Excel files
- Perform kernel density estimation on grain data by sinks and features
- Decompose sink distributions into weighted combination of sources
- Label grains based on learned sources
- Plot kernel densities, sink, and source distributions

```@contents
Depth = 3
```
## How setup the environment
## Temporary method

For v1.1.1, MatrixTensorFactor must be installed separatly first.

`add https://github.com/MPF-Optimization-Laboratory/MatrixTensorFactor.jl.git/#e92375d`
`add https://github.com/njericha/SedimentSourceAnalysis.jl/tree/v1.1.1`

## Prefered Method
1. Run `julia`
2. Add the package with `pkg> add https://github.com/njericha/SedimentSourceAnalysis.jl.git`
(use `julia> ]` to get to the package manager)
3. Resolve dependency issues with `pkg> resolve` and check it works with `pkg> precompile` (if you have an error with dependencies not downloading, try `pkg> instantiate`)
3. Import with `using SedimentSourceAnalysis`

**OR**
### In Browser
1. Go to <https://github.com/njericha/SedimentSourceAnalysis.jl>
2. Click "<> Code" and press "+" to "Create a codespace on main". It make take a few moments to set up.
3. Open the command palett with `Ctrl+Shift+P` (Windows) or `Cmd+Shift+P` (Mac)
4. Enter `>Julia: Start REPL`
5. In the REPL, resolve any dependency issues with `pkg> resolve` and `pkg> instantiate` (use `julia> ]` to get to the package manager). It may take a few minutes to download dependencies.

Run one of the example files by opening the file and pressing the triangular "run" button, or `>Julia: Execute active File in REPL`.

**OR**
### On your own device
1. Clone the repo at <https://github.com/njericha/SedimentSourceAnalysis.jl>
2. Navigate to the root of the repository in a terminal and run `julia`
3. Activate the project with `pkg> activate .` (use `julia> ]` to get to the package manager)
4. resolve any dependency issues with `pkg> resolve`

### Importing the package
Type `julia> using SedimentSourceAnalysis` to load all submodules (currently only `SedimentTools`), or if only one of the modules is desired, type `using SedimentSourceAnalysis.XXX`.

## Examples
`knownsources.jl`: Uses data from Sundel et al where we know the sources of each Grain. Use this to see how well the factorization performs with realistic data.
`unknownsources.jl`: Uses data from Lee et al where we don't have a ground truth. Showcases how the method would be used in practice.
`randomtensor`: Factorizes a random 50x50x50 tensor. See how the factorization performs in theory when a perfect factorization exists.

## Submodules
The main submodule is SedimentTools. SedimentSourceAnalysis also builds off of MatrixTensorFactor.jl.

### MatrixTensorFactor
Defines the main factorization function and related mathematical functions. See the full documentation here [Matrix Tensor Factorization](@ref).

### SedimentTools
Holds various types at the [`Grain`](@ref), and [`Sink`](@ref) level, importing ([`read_raw_data`](@ref)) and processing data ([`make_densities`](@ref)) functions, and additional methods of some [Plots.jl](https://docs.juliaplots.org/stable/) functions for visualization with these custom types.

## Index

```@index
```
