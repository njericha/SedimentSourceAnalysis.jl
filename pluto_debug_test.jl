### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 09c3dbf2-4784-11ee-2a3b-7b08bdb868d5
using Pkg

# ╔═╡ 5468dc35-be18-4a2e-b759-9ab0d14f04c3
Pkg.add("OrderedCollections")

# ╔═╡ 282462fe-060c-44e5-b573-78d8c6f70188
Pkg.add(url="https://github.com/njericha/SedimentSourceAnalysis.jl")

# ╔═╡ c1cc215b-adbc-4ab5-8153-ce1b4aa2b1e3
using SedimentSourceAnalysis

# ╔═╡ 3b65ee86-594a-4bd4-8301-993a4265bd15
using OrderedCollections: OrderedDict

# ╔═╡ d382a339-5540-4381-8ead-ca0600bac34c
begin
	# Import data from excel file
	filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
	sinks = read_raw_data(filename)::Vector{Sink}
end

# ╔═╡ 4a6478a8-ed1e-483a-8d81-4b735df4370f
print(sinks[1][1])

# ╔═╡ Cell order:
# ╠═09c3dbf2-4784-11ee-2a3b-7b08bdb868d5
# ╠═5468dc35-be18-4a2e-b759-9ab0d14f04c3
# ╠═282462fe-060c-44e5-b573-78d8c6f70188
# ╠═c1cc215b-adbc-4ab5-8153-ce1b4aa2b1e3
# ╠═3b65ee86-594a-4bd4-8301-993a4265bd15
# ╠═d382a339-5540-4381-8ead-ca0600bac34c
# ╠═4a6478a8-ed1e-483a-8d81-4b735df4370f
