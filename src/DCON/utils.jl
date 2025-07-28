#This is to implement the original bubble function found in equil/utils.f

# If you are sorting the entire array, just replace the call to bubble with
#
# index = sortperm(key) # for ascending order or
# index = sortperm(key; rev = true) # for descending order

# This will return the indices of the sort range in descending order
function bubble_sorted_indices(key::Vector{Float64}, mmin::Int, mmax::Int)
    inds = mmin:mmax
    sorted_inds = sort(inds; by = i -> key[i], rev = true)
    return collect(sorted_inds)
end

#If you want to mutate an existing index array like your original code:
function bubble!(key::Vector{Float64}, index::Vector{Int}, mmin::Int, mmax::Int)
    sorted_inds = sort(mmin:mmax; by = i -> key[i], rev = true)
    for (j, i) in enumerate(sorted_inds)
        index[mmin + j - 1] = i
    end
end