# Add method to `convert` so I can read strings as Floats
# This should not be used on arbitrary strings and only on numbers surrounded
# by double quotes. This is nessesary because some excel entries have mistakenly
# read the number as a string.
import Base: convert
Base.convert(::Type{Number}, s::String) = parse(T, s)

"""
    data = read_raw_data(filename)

Imports data to a Vector{Sink} type
"""
function read_raw_data(filename)
    # Load the file
    xf = XLSX.readxlsx(filename)

    # Get the list of measurments (each sheet is 1 measurment)
    sheet_names = XLSX.sheetnames(xf)
    isallowed(n) = lowercase(n) ∉ ["source proportions","grain id"]
    filter!(isallowed, sheet_names)

    data = OrderedDict{String, Matrix{Union{Missing,Float64}}}()
    for sheet_name ∈ sheet_names
        @info "extracting $sheet_name..."
        sheet_data = xf[sheet_name][:]
        data[sheet_name] = sheet_data
    end

    data_tensor = cat(collect(values(data)),dims=3)
    n_grains, n_sinks, n_measurements = size(data_tensor)
    list_of_sinks = Vector{Sink}(undef, n_sinks)
    for (j, s) in enumerate(eachslice(data_tensor, dims=(1, 3)))
        vec_of_grains = Vector{<:Grain}(undef, n_grains)
        for (i, g) in enumerate(eachrow(s))
            # Skip grains that contain missing values
            if any(ismissing.(g))
                break
            end
            vec_of_grains[i] = Grain(g::Vector{Float64}, measurment_names=sheet_names)
        end
        list_of_sinks[j] = Sink(vec_of_grains)
    end

    return list_of_sinks
end
