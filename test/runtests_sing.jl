@testset "Test Sing Functions" begin

    # --------------------------
    # Read input matrices from Fortran-compatible .dat files
    # --------------------------
    function read_complex_fortran(fname)
        data = readdlm(fname, Float64)
        m, n = Int(data[1,1]), Int(data[1,2])
        entries = data[2:end, :]
        # Column-major reshape, no transpose
        mat = reshape([Complex(entries[k,1], entries[k,2]) for k in 1:size(entries,1)], m, n)
        return transpose(mat)
    end


    """
        read_solutions_simple(fname) -> Dict{Int, Matrix{Float64}}

    Reads a .dat file with multiple "Solution index" blocks, skipping headers.
    Extracts only the 2nd and 3rd numeric columns and returns a dictionary
    mapping solution index -> (N × 2) matrix.
    """
    function read_solutions_simple(fname::String)
        lines = readlines(fname)
        solutions = Dict{Int, Matrix{Float64}}()

        current_idx = nothing
        col2 = Float64[]
        col3 = Float64[]

        # Helper to store current block when switching solutions
        function flush!()
            if current_idx !== nothing && !isempty(col2)
                solutions[current_idx] = hcat(col2, col3)
                empty!(col2)
                empty!(col3)
            end
        end

        for line in lines
            s = strip(line)
            if isempty(s)
                continue
            end

            # Start of a new solution
            if occursin("Solution index", s)
                flush!()
                current_idx = parse(Int, split(s, "=")[end])
                continue
            end

            # Skip header lines like psifac, q, etc.
            if occursin("=", s)
                continue
            end

            # Split the numeric line by whitespace
            parts = split(s)

            # Expect at least: index col2 col3 ...
            if length(parts) >= 3
                push!(col2, parse(Float64, replace(parts[1+1], 'D' => 'E')))  # column 2
                push!(col3, parse(Float64, replace(parts[1+2], 'D' => 'E')))  # column 3
            end
        end

        # Save the final solution
        flush!()

        print(typeof(solutions), " with ", length(solutions), " solutions read from $fname\n")

        return solutions
    end

    function write_large_matrix(filename, mat::Matrix{ComplexF64})
        m, n = size(mat)
        open(filename, "w") do io
            @printf(io, "%5d %5d\n", m, n)
            for i in 1:m
                for j in 1:n
                    @printf(io, "%16.8e %16.8e\n", real(mat[i,j]), imag(mat[i,j]))
                end
                @printf(io, "----------------------\n")
            end
        end
    end

    # --------------------------
    # Write output in Fortran-compatible format
    # --------------------------
    function write_sing_output(filename, psifac, q, mlow, nn, ud)
        # ud is a 3D array: (mpert, msol, 2), ComplexF64
        mpert, msol, _ = size(ud)

        open(filename, "w") do io
            # header
            @printf(io, "psifac = %16.8f\n", psifac)
            @printf(io, "q = %16.8f\n", q)
            @printf(io, "mlow = %5d\n", mlow)
            @printf(io, "nn = %3d\n", nn)

            # solutions
            for isol in 1:msol
                @printf(io, "Solution index = %3d\n", isol)
                for ipert in 1:mpert
                    xr1 = real(ud[ipert, isol, 1])
                    xi1 = imag(ud[ipert, isol, 1])
                    xr2 = real(ud[ipert, isol, 2])
                    xi2 = imag(ud[ipert, isol, 2])
                    @printf(io, "%5d  %16.8e  %16.8e  %16.8e  %16.8e\n",
                            ipert, xr1, xi1, xr2, xi2)
                end
            end
        end
    end
    
    function read_solutions_3d(fname::String)
        lines = readlines(fname)

        # --- find all solution indices and gather their numeric blocks ---
        blocks = Vector{Vector{Vector{Float64}}}()  # blocks[sol] = rows of 5 columns (as Float64)
        current = Vector{Vector{Float64}}()

        for s in lines
            t = strip(s)
            isempty(t) && continue

            if occursin("Solution index", t)
                if !isempty(current)
                    push!(blocks, current)
                    current = Vector{Vector{Float64}}()
                end
                continue
            end

            if occursin("=", t)  # skip header lines
                continue
            end

            vals = [parse(Float64, replace(x, 'D'=>'E')) for x in split(t)]
            push!(current, vals)
        end
        if !isempty(current)
            push!(blocks, current)
        end

        mpert  = length(blocks[1])           # rows per solution
        nsol   = length(blocks)              # number of solutions
        result = Array{ComplexF64}(undef, mpert, nsol, 2)

        for (j, block) in enumerate(blocks)   # j = solution index
            for (i, row) in enumerate(block)  # i = perturbation index
                # row = [ipert  ud_re ud_im  xss_re xss_im]
                result[i, j, 1] = complex(row[2], row[3])
                result[i, j, 2] = complex(row[4], row[5])
            end
        end
        return result
    end


    #Helper function to extract values like psifac and q from file
    function extract_value(filename::String, item::String)
        for line in eachline(filename)
            parts = split(line, '=')
            if strip(parts[1]) == item
                # remove spaces, carriage returns, and allow D exponents
                val = strip(parts[2]) #replace(strip(parts[2]), r"[dD]" => "E")
                return parse(Float64, val)
            end
        end
        return NaN # fallback if not found
    end

    function copyForSplines(prev_mat, num_range)
        mmat = zeros(ComplexF64, length(num_range), size(prev_mat,1)*size(prev_mat,2))
        #println("size mmat: ", size(mmat))
        #print(size(prev_mat), " ", length(num_range), "\n")
        #print("prev_mat = $(prev_mat)\n")
        for i in 1:length(num_range)
            mmat[i,:] .= vec(prev_mat)
        end
        return mmat
    end

    @testset " Test sing_der " begin
        msol = 32
        mpert = 32
        equil = Equilibrium.setup_equilibrium(joinpath(@__DIR__,"../notebooks/soloviev_ideal/equil.toml"));
        odet = JPEC.DconMod.OdeState(mpert=32, msol=32)
        intr = JPEC.DconMod.DconInternal()
        intr.mpert = msol
        ctrl = JPEC.DconMod.DconControl()

        #=
        # Spline evaluation
        #Calcualte this instead of hardcoding --> import q
        q = extract_value("test_data/sing_der_testing/mat_dat/sing_der_output_normal.dat", "q") #3.26298350
        singfac = mlow .- nn*q .+ collect(0:mpert-1)
        singfac .= 1.0 ./ singfac
        =#

        print(intr.mpert, " ", odet.msol, "\n")

        #TODO: HD5 file eventually? --> look at dcon.jl for examples
        #odet.du_temp = zeros(ComplexF64, intr.mpert, odet.msol, 2)

        psifac_dummy = collect(range(0, 1, 10))
        points = length(psifac_dummy)
        
        amat = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/amat.dat"))
        #println(amat)
        amats = copyForSplines(amat, psifac_dummy)
        write_large_matrix(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/amats_for_splines.dat"), amats)
        #print(amats[])
        bmat = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/bmat.dat"))
        bmats = copyForSplines(bmat, psifac_dummy)
        cmat= read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/cmat.dat"))
        cmats = copyForSplines(cmat, psifac_dummy)
        fmat = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/fmat.dat"))
        fmat .= cholesky(Hermitian(fmat)).L 
        fmats = copyForSplines(fmat, psifac_dummy)
        kmat = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/kmat.dat"))
        kmats = copyForSplines(kmat, psifac_dummy)
        gmat = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/gmat.dat"))
        gmats = copyForSplines(gmat, psifac_dummy)

        umat_p1 = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/umat_p1.dat"))
        umat_p2 = read_complex_fortran(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/umat_p2.dat"))
        odet.psifac = extract_value(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/sing_der_output_normal.dat"), "psifac")

        odet.u[:,:,1] .= umat_p1
        odet.u[:,:,2] .= umat_p2
        #odet.q = extract_value(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/sing_der_output_normal.dat"), "q")#3.26298350
        ctrl.nn = extract_value(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/sing_der_output_normal.dat"), "nn")
        intr.mlow = extract_value(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/sing_der_output_normal.dat"), "mlow")
        #singfac = mlow .- nn*q .+ collect(0:mpert-1)
        #singfac .= 1.0 ./ singfac

        println("nn = $(ctrl.nn), mlow = $(intr.mlow)")
        
        # Initialize structs with relevant data
        ffit = JPEC.DconMod.FourFitVars()
        ffit.amats = SplinesMod.CubicSpline(psifac_dummy, reshape(amats, points, :); bctype="extrap")
        println(size(reshape(amats, points, :)))
        println(typeof(reshape(amats, points, :)))
        write_large_matrix(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/amats_julia.dat"), reshape(amats, points, :))
        println(typeof(ffit.amats))
        ffit.bmats = SplinesMod.CubicSpline(psifac_dummy, reshape(bmats, points, :); bctype="extrap")
        ffit.cmats = SplinesMod.CubicSpline(psifac_dummy, reshape(cmats, points, :); bctype="extrap")
        ffit.fmats = SplinesMod.CubicSpline(psifac_dummy, reshape(fmats, points, :); bctype="extrap")
        ffit.kmats = SplinesMod.CubicSpline(psifac_dummy, reshape(kmats, points, :); bctype="extrap")
        ffit.gmats = SplinesMod.CubicSpline(psifac_dummy, reshape(gmats, points, :); bctype="extrap")

        du = zeros(ComplexF64, intr.mpert, odet.msol, 2)

        # Initialize u with some value, e.g., ones
        # odet.u .= ones(ComplexF64, size(odet.u))
        # Initialize structs with relevant data, can do something like
        # ffit = FourfitVars()
        # amat = load in from test_data folder
        # ffit.amats = SplinesMod.CubicSpline("some psi", reshape(amats, mpsi+1, :); bctype="extrap")
        # if using one psi values causes issues, can make a psi array and just fill it all with amat
        # repeat for other matrices, equil data, other relevant constants, etc.
        
        params = (ctrl, equil, intr, odet, ffit)
        JPEC.DconMod.sing_der!(du, odet.u, params, odet.psifac)

        # load in du from fortran
        # @test du == expected_du (to within some error)
        du_fortran = read_solutions_3d(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/sing_der_output_du.dat")) #sing_der_output_normal.dat holds the ud data; sing_der_output_dat.dat to holds the du data 

        write_sing_output(joinpath(@__DIR__,"test_data/sing_der_testing/mat_dat/sing_der_output_julia.dat"), odet.psifac, odet.q, intr.mlow, ctrl.nn, du)

        println(size(du))
        println(size(du_fortran))
        #println("du = $du")
        #println("du_fortran = $du_fortran")

        @show size(du)
        @show size(du_fortran)
        @show maximum(abs.(du .- du_fortran))
        
        for sol_idx in 1:1#odet.msol
            #du_expected = du_fortran[sol_idx]
            @test isapprox(du[:, sol_idx, 1], du_fortran[:, sol_idx, 1]; rtol=1e-3)
            @test isapprox(du[:, sol_idx, 2], du_fortran[:, sol_idx, 2]; rtol=1e-3)
        end
    end

    @testset " Test sing_find " begin # continue with other functions
    end
end