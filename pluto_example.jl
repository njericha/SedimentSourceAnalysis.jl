### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 05aeff39-5221-4261-9732-a20893f92335
using Pkg

# ╔═╡ 97c93b24-f145-4814-a733-92cca1e9d142
Pkg.add(url="https://github.com/njericha/Sediment-Source-Analysis.jl")

# ╔═╡ f7acb584-bdbc-4fa7-a96a-e5e217454af7
begin
	Pkg.add("XLSX")
	Pkg.add("NamedArrays")
	Pkg.add("Plots")
	Pkg.add("Printf")
	Pkg.add("Random")
	Pkg.add("OrderedCollections")
end

# ╔═╡ da27d62e-84e9-499d-afe1-90ccd61ae5f8
using SedimentAnalysis

# ╔═╡ 0d218304-0758-41d7-8f2f-580c46553769
# note: depending on the environment that Pluto links to (and what packages are already installed there) you may need to also run Pkg.add("package_name") before running this cell.

begin
	using XLSX
	using NamedArrays
	using Plots
	using Printf
	using Random
	#using Statistics: mean
	using OrderedCollections: OrderedDict
end

# ╔═╡ 9eb825f8-1bc7-411e-ac73-9ed77ab69013
begin
	# Import data from excel file
	filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
	sinks = read_raw_data(filename)::Vector{Sink}
end

# ╔═╡ 3ad8323d-f53c-447d-9d1c-3793d707b3c5
begin
	## Look at a grain
	sink1 = sinks[begin]
	grain1 = sink1[begin]
	println("Grain 1 in Sink 1")
	display(grain1)
end

# ╔═╡ Cell order:
# ╠═05aeff39-5221-4261-9732-a20893f92335
# ╠═97c93b24-f145-4814-a733-92cca1e9d142
# ╠═da27d62e-84e9-499d-afe1-90ccd61ae5f8
# ╠═f7acb584-bdbc-4fa7-a96a-e5e217454af7
# ╠═0d218304-0758-41d7-8f2f-580c46553769
# ╠═9eb825f8-1bc7-411e-ac73-9ed77ab69013
# ╠═3ad8323d-f53c-447d-9d1c-3793d707b3c5
