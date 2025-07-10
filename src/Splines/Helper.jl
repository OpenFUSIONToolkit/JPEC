module Helper

export parse_bctype


const BCTYPE_MAP = Dict(
    "natural"      => 1, # unstable
    "periodic"     => 2,
    "extrap" => 3,
    "not-a-knot"   => 4,
    "notaknot" => 4
)

"""
    parse_bctype(bctype)

Internal helper to parse a boundary condition into a validated integer code.

## Arguments:
- `bctype`: The boundary condition as a `String` or `Int`. Valid options are:
    - `"natural"` or `1`
    - `"periodic"` or `2`
    - `"extrap"` or `3`
    - `"not-a-knot"` or `4`

## Returns:
- A validated `Int` code (1-4).
"""
function parse_bctype(bctype::Union{String, Int})
    if bctype isa String

        normalized_bctype = replace(lowercase(bctype), r"[-_\s]" => "")
        if haskey(BCTYPE_MAP, normalized_bctype)
            return BCTYPE_MAP[normalized_bctype]
        else
            error("Invalid string for `bctype`: '$bctype'. Valid options are: $(keys(BCTYPE_MAP))")
        end
    elseif bctype isa Int
        # If it's already an integer, validate it.
        if bctype in values(BCTYPE_MAP)
            return bctype
        else
            error("Invalid integer for `bctype`: $bctype. Valid options are: $(collect(values(BCTYPE_MAP)))")
        end
    end
end


end 