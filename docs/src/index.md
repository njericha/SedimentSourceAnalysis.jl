# Sediment Source Analysis

```@contents
Depth = 3
```
## How setup the environment

### In Browser
1. Go to https://github.com/njericha/Sediment-Source-Analysis.jl
2. Click "<> Code" and press "+" to "Create a codespace on main". It make take a few moments to set up.
3. Open the command palett with `Ctrl+Shift+P` (Windows) or `Cmd+Shift+P` (Mac)
4. Enter `>Julia: Start REPL`
5. In the REPL, resolve any dependency issues with `pkg> resolve` and `pkg> instantiate` (use `julia> ]` to get to the package manager). It may take a few minutes to download dependencies.

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

## Index

```@index
```
