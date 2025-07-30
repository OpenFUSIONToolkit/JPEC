mutable struct FourierSplineType
    handle::Ptr{Cvoid}
    _xs::Vector{Float64}
    _ys::Vector{Float64}
    _fs::Array{Float64, 3}
    mx::Int64
    my::Int64
    mband::Int64
    nqty::Int64
    bctype::Int32
    fit_method::Int32


end

@expose_fields FourierSplineType xs ys fs


function _destroy_fspline(fourier::FourierSplineType)
    if fourier.handle != C_NULL
        ccall((:fspline_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), fourier.handle)
        Core.setfield!(fourier, :handle, C_NULL)
    end
end

function _MakeFourierSpline(mx::Int, my::Int, mband::Int, nqty::Int, bctype::Int32, fit_method::Int32)
    h = Ref{Ptr{Cvoid}}()
    ccall((:fspline_c_create, libspline), Cvoid,
          (Int64, Int64, Int64, Int64, Ref{Ptr{Cvoid}}),
          mx, my, mband, nqty, h)

    handle = h[]
    if handle == C_NULL
        error("Failed to create fspline handle in Fortran library.")
    end

    # Return a partially initialized struct with empty data arrays
    return FourierSplineType(handle,
                             Vector{Float64}(undef, 0),
                             Vector{Float64}(undef, 0),
                             Array{Float64, 3}(undef, 0, 0, 0),
                             mx, my, mband, nqty, bctype, fit_method)
end

function _fspline_setup(xs::Vector{Float64}, ys::Vector{Float64}, fs::Array{Float64, 3}
    , mband::Int, bctype::Int32, fit_method::Int32, fit_flag::Bool)


    mx = length(xs) - 1
    my = length(ys) - 1
    nqty = size(fs, 3)
    fourier = _MakeFourierSpline(mx, my, mband, nqty, Int32(bctype), Int32(fit_method))
    fourier._xs = xs
    fourier._xs = ys
    fourier._fs = fs

    ccall((:fspline_c_setup, libspline), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          fourier.handle, xs, ys, fs)
    if fit_method == 1
        ccall((:fspline_c_fit_1, libspline), Cvoid,
              (Ptr{Cvoid}, Int32, Bool),
              fourier.handle, fourier.bctype, fit_flag)
    elseif fit_method == 2
        ccall((:fspline_c_fit_2, libspline), Cvoid,
              (Ptr{Cvoid}, Int32, Bool),
              fourier.handle, fourier.bctype, fit_flag)
    else
        error("Internal error: Invalid fit_method passed to _fspline_setup.")
    end

    return fourier
end

function fspline_setup(xs::Vector{Float64}, ys::Vector{Float64}, fs::Array{Float64, 3}
    , mband::Int; bctype::Union{String, Int}="not-a-knot", fit_method::Int=1, fit_flag::Bool=true)
    """
    # fspline_setup(xs, ys, fs, mband; bctype="not-a-knot", fit_method=1)

    Creates and fits a function of two variables, f(x, y), to a cubic spline
    in the x-direction and a Fourier series in the y-direction. The y-direction
    is assumed to be periodic.

    ## Arguments:
    - `xs`: Vector of x-coordinates (length `mx`+1).
    - `ys`: Vector of y-coordinates (length `my`+1, periodic direction).
    - `fs`: 3D array of function values with dimensions (`mx`+1, `my`+1, `nqty`).
    - `mband`: Number of Fourier modes (harmonics) to keep, from 0 to `mband`.

    ## Keyword Arguments:
    - `bctype`: Boundary condition type for the cubic spline in `x`.
        - 1: Natural spline (default)
        - 2: Periodic spline
        - 3: Extrapolated spline
        - 4: "Not-a-knot" spline
    - `fit_method`: Algorithm for computing Fourier coefficients.
        - 1: Integration method (for non-uniform `y` grids).
        - 2: Fast Fourier Transform (FFT) method (requires `length(ys)-1` to be a power of 2).

    ## Returns:
    - A `FourierSplineType` object ready for evaluation.
    """
    if length(xs) != size(fs, 1) || length(ys) != size(fs, 2)
        error("Grid vector dimensions must match `fs` array dimensions.")
    end

    if fit_method == 2
        my = length(ys) - 1
        if !ispow2(my)
            error("For `fit_method=2` (FFT), `length(ys)-1` (is $my) must be a power of 2.")
        end
        if 2 * mband > my - 1
            error("For `fit_method=2` (FFT), `2*mband` must not be greater than `(length(ys)-1) - 1`. Got 2*$(mband) > $(my-1).")
        end
    elseif fit_method != 1
        error("Invalid `fit_method`. Choose 1 or 2.")
    end

    local bctype_code::Int = parse_bctype(bctype)

    if bctype_code == 1
        error("Fourier spline  doesn't have natural spline. (bctype = 1/natural)")
    end
    # Call the internal setup function that creates and fits the spline
    fourier = _fspline_setup(xs, ys, fs, mband, Int32(bctype_code), Int32(fit_method),fit_flag)

    # Add a finalizer to ensure the Fortran object is deallocated when the Julia object is garbage collected.
    finalizer(_destroy_fspline, fourier)

    return fourier
end

function _fspline_eval(spl::FourierSplineType, x::Float64, y::Float64, derivs::Int=0)
    if derivs == 0
        f = Vector{Float64}(undef, spl.nqty)
        ccall((:fspline_c_eval, libspline), Cvoid,
            # (handle,      x,       y,      f_out)
            (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}),
            spl.handle, x, y, f)
        return f

    elseif derivs == 1
        f = Vector{Float64}(undef, spl.nqty)
        fx = Vector{Float64}(undef, spl.nqty)
        fy = Vector{Float64}(undef, spl.nqty)
        ccall((:fspline_c_eval_deriv, libspline), Cvoid,
            # (handle,      x,       y,      f_out,      fx_out,      fy_out)
            (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
            spl.handle, x, y, f, fx, fy)
        return f, fx, fy

    elseif derivs == 2
        f = Vector{Float64}(undef, spl.nqty)
        fx = Vector{Float64}(undef, spl.nqty)
        fy = Vector{Float64}(undef, spl.nqty)
        fxx = Vector{Float64}(undef, spl.nqty)
        fxy = Vector{Float64}(undef, spl.nqty)
        fyy = Vector{Float64}(undef, spl.nqty)
        ccall((:fspline_c_eval_deriv2, libspline), Cvoid,
            # (handle,      x,       y,      f_out,      fx_out,      fy_out,      fxx_out,     fxy_out,     fyy_out)
            (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
            spl.handle, x, y, f, fx, fy, fxx, fxy, fyy)
        return f, fx, fy, fxx, fxy, fyy

    else
        error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
    end
end


function _fspline_eval(fspline::FourierSplineType, xs::Vector{Float64}, ys::Vector{Float64}, derivs::Int=0)
	# xs -> Float64 (any length)
	# ys -> Float64 (any length)
	# Returns a matrix of Float64 (length(xs), length(ys), nqty)

	n = length(xs)
	m = length(ys)
	fs = Array{Float64}(undef, n, m, fspline.nqty)
	f = Vector{Float64}(undef, fspline.nqty)
	if derivs > 0
		fsx = Array{Float64}(undef, n, m, fspline.nqty)
		fsy = Array{Float64}(undef, n, m, fspline.nqty)
		fx = Vector{Float64}(undef, fspline.nqty)
		fy = Vector{Float64}(undef, fspline.nqty)
	end
	if derivs > 1
		fsxx = Array{Float64}(undef, n, m, fspline.nqty)
		fsxy = Array{Float64}(undef, n, m, fspline.nqty)
		fsyy = Array{Float64}(undef, n, m, fspline.nqty)
		fxx = Vector{Float64}(undef, fspline.nqty)
		fxy = Vector{Float64}(undef, fspline.nqty)
		fyy = Vector{Float64}(undef, fspline.nqty)
	end
	for i in 1:n
		for j in 1:m
			if derivs == 0
				f = _fspline_eval(fspline, xs[i], ys[j], 0)
				fs[i, j, :] = f
			elseif derivs == 1
				f, fx, fy = _fspline_eval(fspline, xs[i], ys[j], 1)
				fs[i, j, :] = f
				fsx[i, j, :] = fx
				fsy[i, j, :] = fy
			elseif derivs == 2
				f, fx, fy, fxx, fxy, fyy = _fspline_eval(fspline, xs[i], ys[j], 2)
				fs[i, j, :] = f
				fsx[i, j, :] = fx
				fsy[i, j, :] = fy
				fsxx[i, j, :] = fxx
				fsxy[i, j, :] = fxy
				fsyy[i, j, :] = fyy
			else
				error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
			end
		end
	end
	if derivs == 0
		return fs
	elseif derivs == 1
		return fs, fsx, fsy
	elseif derivs == 2
		return fs, fsx, fsy, fsxx, fsxy, fsyy
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
	end
end

function fspline_eval(fourier::FourierSplineType, x, y, derivs::Int=0)
    """
    # fspline_eval(spl, x, y; derivs=0)

    Evaluates a fitted Fourier-Spline at given coordinates.

    ## Arguments:
    - `spl`: A `FourierSplineType` object from `fspline_setup`.
    - `x`: A `Float64` or a `Vector{Float64}` of x-coordinates.
    - `y`: A `Float64` or a `Vector{Float64}` of y-coordinates.

    ## Returns:
    - If `x`, `y` are scalars: A tuple containing the function value(s) and any requested derivatives. If nqty=1, the results are scalars, otherwise they are vectors.
    - If `x`, `y` are vectors: A tuple of 3D arrays for the function values and derivatives on the grid defined by `x` and `y`.
    """
    if derivs < 0 || derivs > 2
        error("Keyword `derivs` must be 0, 1, or 2.")
    end


    results = _fspline_eval(fourier, x, y, derivs)


    if fourier.nqty == 1 && isa(x, Real) && isa(y, Real)
        if isa(results, Tuple)
            return map(vec -> vec[1], results)
        else
            return results[1]
        end
    end

    return results
end