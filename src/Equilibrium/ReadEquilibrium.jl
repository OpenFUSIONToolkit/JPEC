#=
This file contains functions for reading equilibrium files from diferent codes
    that use different formating and collecting the inputs required to form
    a complete PlasmaEquilibrium using either direct or inverse construction
=#


"""
_read_1d_gfile_format(lines_block, num_values)

Internal helper function to parse Fortran-style fixed-width numerical blocks
from a vector of strings.

## Arguments:
- `lines_block`: A `Vector{String}` containing the lines to parse.
- `num_values`: The total number of `Float64` values to read from the block.
## Returns:
- A `Vector{Float64}` containing the parsed values.
"""
function _read_1d_gfile_format(lines_block::Vector{String}, num_values::Int)
    data_str = join(lines_block)
    field_width = 16
    parsed_values = Float64[]
    num_read = 0

    # Ensure the string length is a multiple of the field width for safe processing
    safe_len = (length(data_str) ÷ field_width) * field_width
    for i in 1:field_width:safe_len
        num_read >= num_values && break
        val_str = strip(data_str[i : i + field_width - 1])
        if !isempty(val_str)
            try
                push!(parsed_values, parse(Float64, val_str))
                num_read += 1
            catch e
                @warn "Parsing error for substring: '$val_str'. Error: $e. Skipping."
            end
        end
    end

    if num_read < num_values
        @warn "Expected $num_values values, but only read $num_read."
    end
    return parsed_values
end

"""
    _read_efit(equil_in)

Parses an EFIT g-file, creates initial 1D and 2D splines, and bundles
them into a `DirectRunInput` object.

## Arguments:
- `equil_in`: The `EquilInput` object containing the filename and parameters.
## Returns:
- A `DirectRunInput` object ready for the direct solver.
"""
function read_efit(config::EquilibriumConfig)
    println("--> Processing EFIT g-file: $(config.control.eq_filename)")
    lines = readlines(config.control.eq_filename)

    # --- Parse Header ---
    header1_parts = split(lines[1])
    nw = parse(Int, header1_parts[end-1])
    nh = parse(Int, header1_parts[end])
    println("--> Parsed from header: nw=$nw, nh=$nh")

    header_vals = _read_1d_gfile_format(lines[2:5], 20)
    rdim, zdim, rcentr, rleft, zmid = header_vals[1:5]
    rmaxis, zmaxis, simag, sibry = header_vals[6:9]

    # --- Parse Data Blocks ---
    current_line_idx = 6
    function parse_block(num_pts)
        num_lines = ceil(Int, num_pts / 5)
        block = lines[current_line_idx : current_line_idx + num_lines - 1]
        data = _read_1d_gfile_format(block, num_pts)
        current_line_idx += num_lines
        return data
    end

    fpol_data = parse_block(nw)
    pres_data = parse_block(nw)
    ffprime_data = parse_block(nw)
    pprime_data = parse_block(nw)
    psi_flat_vec = parse_block(nw * nh)
    qprof_data = parse_block(nw)

    psi_rz = reshape(psi_flat_vec, nw, nh)
    println("--> All main data blocks parsed successfully.")

    # --- Create 1D Profile Spline (sq_in) ---
    println("--> Creating 1D profile splines...")
    psi_norm_grid = range(0.0, 1.0, length=nw)
    sq_fs_nodes = hcat(
        abs.(fpol_data),
        max.(pres_data .* mu0, 0.0),
        qprof_data,
        sqrt.(psi_norm_grid)
    )
    # According to Spline_document.txt, bctype=4 is Not-a-Knot
    sq_in = Spl.spline_setup(collect(psi_norm_grid), sq_fs_nodes, bctype=4)
    println("--> 1D Spline fitting complete.")

    # --- Process and Normalize 2D Psi Data ---
    psio_signed = sibry - simag
    psi_proc = (sibry .- psi_rz)
    psio = abs(psio_signed)
    # Ensure psi at the magnetic axis is positive relative to the boundary
    if psio_signed < 0.0
        psi_proc .*= -1.0
    end

    # --- Create 2D Psi Spline (psi_in) ---
    println("--> Creating 2D psi spline...")
    r_grid = range(rleft, rleft + rdim, length=nw)
    z_grid = range(zmid - zdim / 2, zmid + zdim / 2, length=nh)
    rmin, rmax = extrema(r_grid)
    zmin, zmax = extrema(z_grid)

    psi_proc_3d = reshape(psi_proc, (nw, nh, 1))
    psi_in = Spl.bicube_setup(collect(r_grid), collect(z_grid), psi_proc_3d, bctypex=4, bctypey=4)
    println("--> 2D Spline fitting complete.")

    # --- Bundle everything for the solver ---
    return DirectRunInput(config, sq_in, psi_in, rmin, rmax, zmin, zmax, psio)
end

"""
    _read_chease2(config::EquilInput) -> InverseRunInput

Debug version: Reads ASCII CHEASE file and prints checking info at each read.
"""
function read_chease2(config::EquilibriumConfig)
    println("--> Reading CHEASE file: $(config.control.eq_filename)")

    # Robust splitting, also for glued numbers
    split_chease_numbers(str::String) = replace(str, r"([eE][+-]\d+)(?=[\-+]\d+\.)" => s"\1 ") |> split
    parse_chease_floats(str::String) = [parse(Float64, s) for s in split_chease_numbers(str) if !isempty(s)]
    function read_n_floats(io::Any, n::Int; label="unknown")
        vals = Float64[]
        read_lines = 0
        while length(vals) < n && !eof(io)
            line = readline(io)
            read_lines += 1
            newvals = parse_chease_floats(line)
            if !isempty(newvals)
                append!(vals, newvals)
            end
        end
        println("    [DEBUG] $label: read $n floats, actually got $(length(vals)); lines read = $read_lines")
        length(vals) < n && error("[DEBUG] Unexpected EOF while reading $label — needed $n, got $(length(vals))")
        return vals[1:n]
    end
    function read_2d(io::Any, rows::Int, cols::Int; label="unknown")
        arr = Array{Float64}(undef, rows, cols)
        for j in 1:cols
            arr[:,j] .= read_n_floats(io, rows; label="$label col $j")
        end
        println("    [DEBUG] $label: finished 2D read of ($rows x $cols) values")
        arr
    end

    open(config.control.eq_filename, "r") do io
        # Header and axx
        header_line = readline(io)
        header = split(strip(header_line))
        ntnova, npsi1, nsym = parse.(Int, header)
        println("--> Header: ntnova=$ntnova, npsi1=$npsi1, nsym=$nsym")
        axxline = readline(io)
        println("    [DEBUG] axx raw: $axxline")
        axx = parse_chease_floats(axxline)
        println("    [DEBUG] axx parsed: $axx")
        ma = npsi1 - 1
        nrc = ntnova + 3

        # 1D arrays
        zcpr   = read_n_floats(io, ma;      label="zcpr  (ma = $ma)")
        zcppr  = read_n_floats(io, npsi1;   label="zcppr (npsi1 = $npsi1)")
        zq     = read_n_floats(io, npsi1;   label="zq")
        zdq    = read_n_floats(io, npsi1;   label="zdq")
        ztmf   = read_n_floats(io, npsi1;   label="ztmf")
        ztp    = read_n_floats(io, npsi1;   label="ztp")
        zfb    = read_n_floats(io, npsi1;   label="zfb")
        zfbp   = read_n_floats(io, npsi1;   label="zfbp")
        zpsi   = read_n_floats(io, npsi1;   label="zpsi")
        zpsim  = read_n_floats(io, ma;      label="zpsim")

        # 2D arrays
        println("    [DEBUG] About to read 2D: zrcp nrc=$nrc npsi1=$npsi1 (should be $((nrc)*(npsi1)) numbers)")
        zrcp  = read_2d(io, nrc, npsi1;     label="zrcp")
        println("    [DEBUG] About to read 2D: zzcp")
        zzcp  = read_2d(io, nrc, npsi1;     label="zzcp")
        println("    [DEBUG] About to read 2D: zjacm")
        zjacm = read_2d(io, nrc, npsi1;     label="zjacm")
        println("    [DEBUG] About to read 2D: zjac")
        zjac  = read_2d(io, nrc, npsi1;     label="zjac")

        ro = zrcp[1, 1]
        zo = zzcp[1, 1]
        psio = abs(zpsi[end] - zpsi[1])

        println("--> Finished reading CHEASE equilibrium.")
        println("    Magnetic axis at (ro=$ro, zo=$zo), psio=$psio")
        return nothing # replace with construction of your InverseRunInput, if desired
    end
end

