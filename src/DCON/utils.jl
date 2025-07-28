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