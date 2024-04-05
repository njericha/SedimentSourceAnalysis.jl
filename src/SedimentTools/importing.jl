# Add method to `convert` so I can read strings as Floats
# This should not be used on arbitrary strings and only on numbers surrounded
# by double quotes. This is nessesary because some excel entries have mistakenly
# read the number as a string.
#import Base: convert
#Base.convert(::Type{Number}, s::String) = parse(T, s)

"""
    read_raw_data(filename)
    read_raw_data(filename; skip_sheets)

Imports excel data to a `Vector{Sink}`.

Excel file must have one element per page where different columns correspond to different
sinks. Each sink can have a different number of grains (length of the column), but a sink
must have a consistant length across different measurements (sheets).

Optionaly provide a collection `skip_sheets` to blacklist sheet names from the excel file.
"""
function read_raw_data(filename; skip_sheets=Set(["source proportions","grain id"]))
    # Load the file
    xf = readxlsx(filename)

    # Get the list of measurements (each sheet is 1 measurement)
    sheet_names = sheetnames(xf)
    isallowed(n) = lowercase(n) ∉ skip_sheets
    filter!(isallowed, sheet_names)

    # Extract data to a dictionary Dict("element 1" => sheet_data, ...)
    data = OrderedDict{String, Matrix{Union{Missing,Float64}}}()
    for sheet_name ∈ sheet_names
        @info "extracting $sheet_name..."
        sheet_data = xf[sheet_name][:]
        data[sheet_name] = sheet_data
    end

    @assert allequal([size(sheet) for sheet in values(data)])

    # Collect all sheets into a single 3-tensor (Array{T, 3})
    data_tensor = cat(values(data)...,dims=3)
    n_grains, n_sinks, n_measurements = size(data_tensor)

    # Construct a vector of sinks, each sink is a vector of grains
    vec_of_sinks = Vector{Sink}(undef, n_sinks)
    for (j, s) in enumerate(eachslice(data_tensor, dims=2))

        vec_of_grains = Grain[] #Vector{Grain}(Grain(Float64[],String[]), n_grains)
        for (i, g) in enumerate(eachrow(s))

            # Skip grains that contain missing values
            if any(ismissing.(g))
                vec_of_grains = vec_of_grains[1:(i-1)]
                break
            else
                g = Vector{Float64}(g)
            end

            #vec_of_grains[i] = Grain(g, sheet_names) # measurement_names = sheet_names
            push!(vec_of_grains, Grain(g, sheet_names))
        end

        vec_of_sinks[j] = Sink(vec_of_grains)
    end

    return vec_of_sinks
end
