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

function dump_matrix(filename, mat)
    open(filename, "w") do io
        header = ["i", "j", "Re(val)", "Im(val)"]
        println(io, join(header, "\t"))
        ni, nj = size(mat)
        for i in 1:ni
            for j in 1:nj
                println(io, join([i, j, real(mat[i, j]), imag(mat[i, j])], "\t"))
            end
        end
    end
end

function dump_u_matrix(filename, mat)
    open(filename, "w") do io
        header = ["i", "j", "Re(val1)", "Im(val1)", "Re(val2)", "Im(val2)"]
        println(io, join(header, "\t"))
        ni, nj = size(mat)
        for i in 1:ni
            for j in 1:nj
                println(io, join([i, j, real(mat[i, j, 1]), imag(mat[i, j, 1]), real(mat[i, j, 2]), imag(mat[i, j, 2])], "\t"))
            end
        end
    end
end

function dump_matrix_3D(filename, mat)
    open(filename, "w") do io
        header = ["psi", "i", "j", "Re(val)", "Im(val)"]
        println(io, join(header, "\t"))
        mpsi, mpert, _ = size(mat)
        for ipsi in 1:mpsi     
            for j in 1:mpert
                for i in 1:mpert
                    println(io, join([ipsi, i, j, real(mat[ipsi, i, j]), imag(mat[ipsi, i, j])], "\t"))
                end
            end
        end
    end
end

function load_u_matrix(filename)
    lines = readlines(filename)
    data = [parse.(Float64, split(l)) for l in lines[2:end]]
    i_vals = Int.(getindex.(data, 1))
    j_vals = Int.(getindex.(data, 2))
    ncols = (length(data[1]) - 2) ÷ 2
    imax = maximum(i_vals)
    jmax = maximum(j_vals)
    mat = zeros(ComplexF64, imax, jmax, ncols)
    for row in data
        i = Int(row[1])
        j = Int(row[2])
        for k in 1:ncols
            re = Float64(row[2*k+1])
            im = Float64(row[2*k+2])
            mat[i, j, k] = complex(re, im)
        end
    end
    return mat
end