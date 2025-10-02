# This is to implement the original bubble function found in equil/utils.f
# If you are sorting the entire array, just replace the call to bubble with
#
# index = sortperm(key) # for ascending order or
# index = sortperm(key; rev = true) # for descending order

# This will return the indices of the sort range in descending order
function sortperm_subrange(key::Vector{Float64}, mrange::UnitRange{Int}; descend = true )
    sorted_inds = sort(mrange; by = i -> key[i], rev = descend)
    return collect(sorted_inds)
end

#If you want to mutate an existing index array like your original code:
function sortperm_subrange!(key::Vector{Float64}, index::Vector{Int}, mrange::UnitRange{Int}; descend = true )
    sorted_inds = sort(mrange; by = i -> key[i], rev = descend)
    for (j, i) in enumerate(sorted_inds)
        index[first(mrange) + j - 1] = i
    end
end

function load_u_matrix(filename)
    lines = readlines(filename)
    data = [parse.(Float64, split(l)) for l in lines[2:end]]
    i_vals = Int.(getindex.(data, 1))
    j_vals = Int.(getindex.(data, 2))
    k_vals = Int.(getindex.(data, 3))
    imax = maximum(i_vals)
    jmax = maximum(j_vals)
    kmax = maximum(k_vals)
    mat = zeros(ComplexF64, imax, jmax, kmax)
    for row in data
        i = Int(row[1])
        j = Int(row[2])
        k = Int(row[3])
        re = Float64(row[4])
        im = Float64(row[5])
        mat[i, j, k] = complex(re, im)
    end
    return mat
end

"""
    init_files!(df::DconFileNames; append=false, dir=nothing)

Open all files in the struct and store the handles in `df.handles`.
If `append=true`, open in append mode; otherwise overwrite.
If `dir` is provided, files will be created in that directory.
"""
function init_files!(df::DconFileNames; append=false, dir::Union{Nothing,String}=nothing)
    mode = append ? "a" : "w"

    for field in fieldnames(DconFileNames)
        if field == :handles
            continue
        end
        filename = getfield(df, field)

        # if dir is given, join it with the filename
        fullpath = isnothing(dir) ? filename : joinpath(dir, filename)

        # make sure the directory exists
        if !isnothing(dir)
            mkpath(dir)
        end
        io = open(fullpath, mode)
        df.handles[String(field)] = io
    end
    return df
end

"""
    write_to!(df::DconFileNames, field::Symbol, data::AbstractString)

Write data to the file associated with `field` in `df`.

- If all `args` are strings, they are joined with spaces and written
  as a line of text (terminated with a newline).
- Otherwise, `args` are written in sequence as raw binary values,
  matching Fortran-style `WRITE` behavior.
"""
function write_to!(df::DconFileNames, field::Symbol, args...)
    io = df.handles[String(field)]

    if all(x -> x isa AbstractString, args)
        # if all arguments are strings, treat as text line
        println(io, join(args, " "))
    else
        # otherwise, dump binary representations
        write(io, args...)
    end
    flush(io)
end

"""
    close_files!(df::DconFileNames)

Close all open file handles.
"""
function close_files!(df::DconFileNames)
    for io in values(df.handles)
        close(io)
    end
    empty!(df.handles)
end