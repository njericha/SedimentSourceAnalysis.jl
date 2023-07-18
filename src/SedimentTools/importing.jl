# Add method to `convert` so I can read strings as Floats
# This should not be used on arbitrary strings and only on numbers surrounded
# by double quotes. This is nessesary because some excel entries have mistakenly
# read the number as a string.
import Base: convert
Base.convert(::Type{Number}, s::String) = parse(T, s)

"""
    grain_data, rock_index = read_xlsx_data(elements)

Read the excel files for the elements given and return a matrix grain_data where
- each row is data for the same grain
- each column is data for the same element
and rock_index which gives the starting row index for each rock sample.
"""
function read_xlsx_data(generic_path, elements)
    # Initialize empty vectors
    ms,ns,missing_indexs,number_of_grains_per_rocks,rock_indexs,element_cols=([] for _ ∈ 1:6)

    # Extract element data from each file into a long column
    for element ∈ elements
        # Read data
        element_data = XLSX.readdata(generic_path(element),"Sheet1","A2:X61")
        element_data = element_data[:,begin:2:end] # remove every other column to clear "999" columns
        m,n = size(element_data)

        # Find length of data in each column
        missing_index = [findfirst(ismissing.(col)) for col ∈ eachcol(element_data)]
        number_of_grains_per_rock = [(isnothing(i) ? m : i-1) for i ∈ missing_index]
        calc_rock_index(i::Integer) = 1 + sum(number_of_grains_per_rock[1:i-1])
        rock_index = zip(calc_rock_index.(1:n), calc_rock_index.(1:n) .+ number_of_grains_per_rock .- 1) |> collect

        # Reshape data into a long column
        element_col = vcat((col[begin:i] for (col, i) ∈ zip(eachcol(element_data), number_of_grains_per_rock))...)

        # Push all info into a vector
        push!(ms, m)
        push!(ns, n)
        push!(missing_indexs, missing_index)
        push!(number_of_grains_per_rocks, number_of_grains_per_rock)
        push!(rock_indexs, rock_index)
        push!(element_cols, element_col)
    end

    # Check all sizes and grains are equal for each element
    for v ∈ (ms, ns, missing_indexs, number_of_grains_per_rocks, rock_indexs)
        @assert allequal(v) "The list $v does not contain identical elements"
    end

    # Stack columns into a matrix
    grain_data = hcat(element_cols...)

    # Check for negative data and flip them to positive if so
    if !all(grain_data .>= 0)
        negative_indexes = findall(grain_data .< 0)
        @warn "Negative data was found in $negative_indexes, flipping values to be positive."
        grain_data[negative_indexes] .*= -1
    end

    return grain_data, rock_indexs[1]
end

function import_grain_data(filepath, elements)
    # Initialize empty vectors
    ms,ns,missing_indexs,number_of_grains_per_rocks,rock_indexs,element_cols=([] for _ ∈ 1:6)

    # Extract element data from each file into a long column
    for element ∈ elements
        # Read data
        element_data = XLSX.readdata(filepath,element,"A1:T75")
        m,n = size(element_data)

        # Find length of data in each column
        missing_index = [findfirst(ismissing.(col)) for col ∈ eachcol(element_data)]
        number_of_grains_per_rock = [(isnothing(i) ? m : i-1) for i ∈ missing_index]
        calc_rock_index(i::Integer) = 1 + sum(number_of_grains_per_rock[1:i-1])
        rock_index = zip(calc_rock_index.(1:n), calc_rock_index.(1:n) .+ number_of_grains_per_rock .- 1) |> collect

        # Reshape data into a long column
        element_col = vcat((col[begin:i] for (col, i) ∈ zip(eachcol(element_data), number_of_grains_per_rock))...)

        # Push all info into a vector
        push!(ms, m)
        push!(ns, n)
        push!(missing_indexs, missing_index)
        push!(number_of_grains_per_rocks, number_of_grains_per_rock)
        push!(rock_indexs, rock_index)
        push!(element_cols, element_col)
    end

    # Check all sizes and grains are equal for each element
    for v ∈ (ms, ns, missing_indexs, number_of_grains_per_rocks, rock_indexs)
        @assert allequal(v) "The list $v does not contain identical elements"
    end

    # Stack columns into a matrix
    grain_data = hcat(element_cols...)

    # Check for negative data and flip them to positive if so
    if !all(grain_data .>= 0)
        negative_indexes = findall(grain_data .< 0)
        @warn "Negative data was found in $negative_indexes, flipping values to be positive."
        grain_data[negative_indexes] .*= -1
    end

    return grain_data, rock_indexs[1]
end

###

function read_distribution_data(filepath)
    xf = XLSX.readxlsx(filepath)
    sheet_names = XLSX.sheetnames(xf)
    distributions = DistributionData[]
    for sheet_name ∈ sheet_names
        @info "extracting $sheet_name..."
        sheet = XLSX.readtable(filepath, sheet_name)
        sheet_data = hcat(sheet.data...)
        #sheet_data = xf[sheet_name][:]
        scale, data = sheet_data[:,1], sheet_data[:,2:end]

        # Remove columns if they are entirely missing
        if any(ismissing.(data))
            data = hcat(filter(v->(!all(ismissing.(v))), eachcol(data))...)
        end

        data = convert(Matrix{Float64},data)
        scale = convert(Vector{Float64},scale)

        #check for negatives and NaN's
        if !all(data .>= 0)
            negative_indexes = findall(data .< 0)
            n_negative_idx = length(negative_indexes)
            @warn "Negative data was found in $n_negative_idx location(s), flipping values to be positive."
            data[negative_indexes] .*= -1
        end

        if any(data .=== NaN)
            nan_indexes = findall(data .=== NaN)
            n_nan_idx = length(negative_indexes)
            @warn "NaN data was found in $n_nan_idx location(s). Setting them to 0."
            data[nan_indexes] .= 0
        end

        distribution_data = DistributionData(sheet_name, scale, data)
        push!(distributions, distribution_data) #TODO optimize collecting data into a vector
    end
    return Distributions(distributions)
end

"""
    data = read_raw_data(filename)

imports data to an OrderedDict
ex. data == OrderedDict("ages" => [1 2 3 ; 4 5 6], ...)
"""
function read_raw_data(filename)
    xf = XLSX.readxlsx(filename)
    sheet_names = XLSX.sheetnames(xf)
    isallowed(n) = lowercase(n) ∉ ["source proportions","grain id"]
    measurments = filter(isallowed, sheet_names)
    data = OrderedDict{String, Matrix{Union{Missing,Float64}}}()
    for sheet_name ∈ measurments
        @info "extracting $sheet_name..."
        sheet_data = xf[sheet_name][:]
        data[sheet_name]= sheet_data
    end
    return data
end
