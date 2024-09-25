# This script is used to check if the measurements within a source are independent
# This checks the validity of representing the grain distributions as a product
# distribution.

using XLSX
using NamedArrays
using Plots
using StatsPlots

gr(size = (1000, 1000)) # need large figure size to fit correlation plot

using MatrixTensorFactor
using SedimentSourceAnalysis

# Load source data
filename = "./data/sundell2022/3Sources from Sundell et al 2022.xlsx"
#filename = "./data/sundell2022/20sinks from 3Sources from Sundell et al 2022.xlsx"
#filename = "./data/lee2021/Lee et al 2021 All Measurements.xlsx"
sources = read_raw_data(filename)::Vector{Source};

# Combine all grains into a 7×614 matrix
# Column vectors are the grains, and each row is a different measurement
M = hcat((sources...)...)
display(corrplot(M', labels=getmeasurements(sources[1])))

# Mathematical Examples
using Random
n = 5000

# Independent variables
u = randn(n)
v = randexp(n)
M = hcat(u,v)
display(corrplot(M))

u = randn(n)
v = randn(n)
M = hcat(u,v)
display(corrplot(M))

# Dependent variables
u = randn(n)
v = @. exp(-u) + $(randn(n))
M = hcat(u,v)
display(corrplot(M)) # u and v are correlated and dependent

u = randn(n)
v = @. exp(-abs(u)) + 0.2*$(randn(n))
M = hcat(u,v)
display(corrplot(M)) # u and v are NOT correlated but ARE dependent
