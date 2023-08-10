### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 90219b4a-b620-488c-bdb6-d4533721a5da
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

    import SedimentAnalysis
	#import SedimentTools
end

# ╔═╡ 9ce69737-69a1-4577-9a96-359b3d580fbf
begin
	using XLSX: readxlsx, sheetnames
	using NamedArrays: NamedArray, NamedMatrix, NamedVector, setnames!, setdimnames!, dimnames
	using Plots
	using OrderedCollections: OrderedDict
	#using Base: AbstractVecOrTuple
end

# ╔═╡ 9a468615-3ce4-4b33-a9e8-5e9f89f01e7b
include("./src/SedimentTools/importing.jl")

# ╔═╡ ec3fd2dc-a07e-4615-8f55-aaf68e8f0d73
include("./src/SedimentTools/structs.jl")

# ╔═╡ 464e30a2-38bb-4b73-896a-4e7a162a86ce
# ╠═╡ disabled = true
#=╠═╡
#import Pkg
  ╠═╡ =#

# ╔═╡ f8b90b60-ae77-4041-9a41-a529e54fbfda
# ╠═╡ disabled = true
#=╠═╡
include("./src/SedimentTools/SedimentTools.jl")
#/Users/naomi/Desktop/NTD/Sediment-Source-Analysis.jl/src/SedimentTools
#/Users/naomi/Desktop/NTD/Sediment-Source-Analysis.jl/src/SedimentTools.jl
  ╠═╡ =#

# ╔═╡ 3ad8323d-f53c-447d-9d1c-3793d707b3c5


# ╔═╡ d18c3684-dd84-4b53-af7a-5c4b87cea70b
# ╠═╡ disabled = true
#=╠═╡
#using SedimentTools
#Pkg.add("SedimentTools")
  ╠═╡ =#

# ╔═╡ 1b8c5969-65fe-4175-88b7-b76d80063a23
begin
	# Import data from excel file
	filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
	sinks = read_raw_data(filename)::Vector{Sink}
end

# ╔═╡ cb573e26-5068-4105-a307-a50fcf1b20c3
# ╠═╡ disabled = true
#=╠═╡
#include("./src/SedimentAnalysis.jl")
  ╠═╡ =#

# ╔═╡ 4a9d6b6a-36d8-11ee-2c7d-fbc5a9b6bd29
# ╠═╡ disabled = true
#=╠═╡

#begin
#	using XLSX
#	using NamedArrays
#	using Plots
#end

  ╠═╡ =#

# ╔═╡ cec5b85d-1a48-4c8a-9cce-d4b4e2bd702d
# ╠═╡ disabled = true
#=╠═╡
#begin
#	Pkg.add("XLSX")
#	Pkg.add("NamedArrays")
#	Pkg.add("Plots")
#end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═464e30a2-38bb-4b73-896a-4e7a162a86ce
# ╠═90219b4a-b620-488c-bdb6-d4533721a5da
# ╠═f8b90b60-ae77-4041-9a41-a529e54fbfda
# ╠═9ce69737-69a1-4577-9a96-359b3d580fbf
# ╠═3ad8323d-f53c-447d-9d1c-3793d707b3c5
# ╠═9a468615-3ce4-4b33-a9e8-5e9f89f01e7b
# ╠═ec3fd2dc-a07e-4615-8f55-aaf68e8f0d73
# ╠═d18c3684-dd84-4b53-af7a-5c4b87cea70b
# ╠═1b8c5969-65fe-4175-88b7-b76d80063a23
# ╠═cb573e26-5068-4105-a307-a50fcf1b20c3
# ╠═4a9d6b6a-36d8-11ee-2c7d-fbc5a9b6bd29
# ╠═cec5b85d-1a48-4c8a-9cce-d4b4e2bd702d
