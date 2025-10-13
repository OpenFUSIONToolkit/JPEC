# This is to implement the original bubble function found in equil/utils.f
# If you are sorting the entire array, just replace the call to bubble with
#
# index = sortperm(key) # for ascending order or
# index = sortperm(key; rev = true) # for descending order

# This will return the indices of the sort range in descending order
function sortperm_subrange(key::Vector{Float64}, mrange::UnitRange{Int}; descend=true)
    sorted_inds = sort(mrange; by=i -> key[i], rev=descend)
    return collect(sorted_inds)
end

#If you want to mutate an existing index array like your original code:
function sortperm_subrange!(key::Vector{Float64}, index::Vector{Int}, mrange::UnitRange{Int}; descend=true)
    sorted_inds = sort(mrange; by=i -> key[i], rev=descend)
    for (j, i) in enumerate(sorted_inds)
        index[first(mrange)+j-1] = i
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
    init_files(out::DconOutput, path::String)

Open all requested output files into `path` and store handles.
"""
function init_files(out::DconOutput, path::String)
    for pname in propertynames(out)
        if startswith(String(pname), "write_") && getproperty(out, pname)
            base_str = String(pname)[7:end]          # remove write_ prefix
            base = Symbol(base_str)                       # key for handles
            fname_field = Symbol("fname_" * base_str) # add fname_ prefix
            fname = getproperty(out, fname_field)

            full_path = joinpath(path, fname)

            if endswith(fname, ".out") || endswith(fname, ".txt")
                out.handles[base] = open(full_path, "w")

            elseif endswith(fname, ".h5") && !haskey(out.handles, base)
                fid = h5open(full_path, "w")
                out.handles[base] = fid

            else
                error("Unknown file extension for $fname, need to add to init_files!")
            end
        end
    end
    return out
end

"""
    write_output(out::DconOutput, key::Symbol, data; dsetname=nothing, slice=:)

Writes `data` to the file corresponding to `key`.

  - For text files: appends a line.
  - For HDF5 files: writes into `dsetname` at `slice`.
"""
function write_output(out::DconOutput, key::Symbol, data; dsetname=nothing, slice=:)
    handle = out.handles[key]

    if handle isa IOStream
        println(handle, data)

    elseif handle isa HDF5.File
        if dsetname === nothing
            error("Must specify dsetname when writing to HDF5 file $key")
        end
        if haskey(handle, dsetname)
            handle[dsetname][slice...] = data
        else
            # create dataset if not present
            handle[dsetname] = data
        end

    else
        error("Unsupported handle type for key $key")
    end
end

"""
    close_files(out::DconOutput)

Closes all open files.

# TODO: is there a way to make the code run this if it errors out anywhere?

# In existing state, files will remain open if there is an error, so you have

# to manually exit and reopen Julia.
"""
function close_files(out::DconOutput)
    for (key, handle) in out.handles
        close(handle)
    end
    return empty!(out.handles)
end

"""
    chebyshev_nodes(a::Float64, b::Float64, N::Int)
Generates `N` Chebyshev-Lobatto nodes in the interval `[a, b]`
in ascending order.
"""
# TODO: this is no longer used, but might be useful code? Leaving for now, but likely can be removed
function chebyshev_nodes(a::Float64, b::Float64, N::Int)
    j = 0:N-1
    nodes = (a + b) / 2 .+ (b - a) / 2 .* cos.(π * j ./ (N - 1))
    return reverse(nodes)
end