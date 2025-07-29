module CubicSpline

const libdir = joinpath(@__DIR__, "..", "..", "deps")
const libspline = joinpath(libdir, "libspline")

using ..Helper: parse_bctype, ReadOnlyArray, @expose_fields

export spline_setup, spline_eval, spline_eval!, spline_integrate!

abstract type CubicSplineType end

mutable struct RealSplineType <: CubicSplineType
	handle::Ptr{Cvoid}
	_xs::Vector{Float64} #x-coordinate (size mx)
	_fs::Matrix{Float64} #spline values and their derivatives (size mx*nqty)
	mx::Int64 #number of coordinate values in _xs
	nqty::Int64 #number of quantities in _fs
    bctype::Int32  # Boundary condition type

	_fsi::Matrix{Float64} # To store integrals at gridpoint
	_fs1::Matrix{Float64} # To store 1-deriv at gridpoint

end

mutable struct ComplexSplineType <: CubicSplineType
	handle::Ptr{Cvoid}
	_xs::Vector{Float64} #x-coordinate (size mx)
	_fs::Matrix{ComplexF64} #spline values and their derivatives (size mx*nqty)
	mx::Int64 #number of coordinate values in _xs
	nqty::Int64 #number of quantities in _fs
    bctype::Int32  # Boundary condition type

	_fsi::Matrix{ComplexF64} # To store integrals at gridpoint
	_fs1::Matrix{ComplexF64} # To store 1-deriv at gridpoint

end

@expose_fields RealSplineType xs fs fsi fs1
@expose_fields ComplexSplineType xs fs fsi fs1


function _destroy_spline(spline::CubicSplineType)
    if spline.handle != C_NULL
        ccall((:spline_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), spline.handle)
		Core.setfield!(spline, :handle, C_NULL)
    end
end

function _destroy_spline(spline::ComplexSplineType)
    if spline.handle != C_NULL
        ccall((:cspline_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), spline.handle)
		Core.setfield!(spline, :handle, C_NULL)
    end
end


function _MakeSpline(mx::Int64, nqty::Int64, bctype::Int32)
	h = Ref{Ptr{Cvoid}}()
	ccall((:spline_c_create, libspline), Cvoid,
		(Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)
	

	fsi = Matrix{Float64}(undef, 0, 0)
	fs1 = Matrix{Float64}(undef, 0, 0)
	
	return RealSplineType(h[], Vector{Float64}(undef, 0), Matrix{Float64}(undef, 0, 0)
	, mx, nqty, bctype, fsi, fs1)
end

function _MakeCSpline(mx::Int64, nqty::Int64, bctype::Int32)
	h = Ref{Ptr{Cvoid}}()
	ccall((:cspline_c_create, libspline), Cvoid,
		(Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)

	fsi = Matrix{ComplexF64}(undef, 0, 0)
	fs1 = Matrix{ComplexF64}(undef, 0, 0)

	return ComplexSplineType(h[], Vector{Float64}(undef, 0), Matrix{ComplexF64}(undef, 0, 0),
	 mx, nqty, bctype, fsi, fs1)
end





function _spline_setup(xs::Vector{Float64}, fs::Vector{Float64}, bctype::Int32)
	# xs -> Float64 (mx)
	# fs -> Float64 (mx, nqty)
	if length(xs) != length(fs)
		error("Length of xs must match length of fs")
	end
	mx = length(xs)-1
	nqty = 1  # Default to 1 quantity if not specified
	spline = _MakeSpline(mx, nqty, Int32(bctype))
	spline._xs = xs
	# Convert fs to a matrix with one column
	fs_matrix = reshape(fs, mx+1, nqty) #considering definition, we need mx+1
	spline._fs = fs_matrix
	spline._fs1 = Matrix{Float64}(undef, mx+1, nqty)

	ccall((:spline_c_setup, libspline), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}),
		spline.handle, xs, fs_matrix)

	ccall((:spline_c_fit, libspline), Cvoid, 
        (Ptr{Cvoid}, Int32, Ptr{Float64}), spline.handle, spline.bctype, spline._fs1)


	return spline
end

function _spline_setup(xs::Vector{Float64}, fs::Matrix{Float64}, bctype::Int32)
	# xs -> Float64 (mx)
	# fs -> Float64 (mx, nqty)
	if length(xs) != size(fs, 1)
		error("Length of xs must match number of rows in fs")
	end
	mx = length(xs)-1
	nqty = size(fs, 2)
	spline = _MakeSpline(mx, nqty, Int32(bctype))
	spline._xs = xs
	spline._fs = fs
	spline._fs1 = Matrix{Float64}(undef, mx+1, nqty)

	ccall((:spline_c_setup, libspline), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}),
		spline.handle, xs, fs)

	ccall((:spline_c_fit, libspline), Cvoid, 
        (Ptr{Cvoid}, Int32, Ptr{Float64}), spline.handle, spline.bctype, spline._fs1)
	return spline
end

function _spline_setup(xs::Vector{Float64}, fs::Vector{ComplexF64}, bctype::Int32)
	# xs -> Float64 (mx)
	# fs -> ComplexF64 (mx, nqty)
	if length(xs) != length(fs)
		error("Length of xs must match length of fs")
	end
	mx = length(xs)-1
	nqty = 1  # Default to 1 quantity if not specified
	spline = _MakeCSpline(mx, nqty, Int32(bctype))
	spline._xs = xs
	# Convert fs to a matrix with one column
	fs_matrix = reshape(fs, mx+1, nqty)
	spline._fs = fs_matrix
	spline._fs1 = Matrix{ComplexF64}(undef, mx+1, nqty)

	ccall((:cspline_c_setup, libspline), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{ComplexF64}),
		spline.handle, xs, fs_matrix)

    ccall((:cspline_c_fit, libspline), Cvoid,
        (Ptr{Cvoid}, Int32, Ptr{ComplexF64}), spline.handle, spline.bctype, spline._fs1)
	return spline
end

function _spline_setup(xs::Vector{Float64}, fs::Matrix{ComplexF64}, bctype::Int32)
	# xs -> Float64 (mx)
	# fs -> ComplexF64 (mx, nqty)
	if length(xs) != size(fs, 1)
		error("Length of xs must match number of rows in fs")
	end
	mx = length(xs)-1
	nqty = size(fs, 2)
	spline = _MakeCSpline(mx, nqty, Int32(bctype))
	spline._xs = xs
	spline._fs = fs
	spline._fs1 = Matrix{ComplexF64}(undef, mx+1, nqty)


	ccall((:cspline_c_setup, libspline), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{ComplexF64}),
		spline.handle, xs, fs)
    ccall((:cspline_c_fit, libspline), Cvoid,
		(Ptr{Cvoid}, Int32, Ptr{ComplexF64}), spline.handle, spline.bctype, spline._fs1)
	return spline
end

function spline_setup(xs, fs; bctype::Union{String, Int}="not-a-knot")
    """
    # spline_setup(xs, fs, bctype="not-a-knot")
    ## Arguments:
    - `xs`: A vector of Float64 values representing the x-coordinates.
    - `fs`: A vector or matrix of Float64/ComplexF64 values representing the function values at the x-coordinates.
    ## Keyword Arguments:
    - `bctype`: Boundary condition type for the cubic spline in `x`.
        - 1: Natural spline (default)
        - 2: Periodic spline
        - 3: Extrapolated spline
        - 4: "Not-a-knot" spline
    ## Returns:
    - A `Spline` object containing the spline handle, x-coordinates, function values,
      number of x-coordinates, number of quantities, and index of x position in the spline.
    """
	if !isa(xs, Vector{Float64})
		error("xs must be a vector of Float64")
	end
	if !isa(fs, Vector{Float64}) && !isa(fs, Matrix{Float64}) && !isa(fs, Vector{ComplexF64}) && !isa(fs, Matrix{ComplexF64})
		error("fs must be a vector or matrix of Float64 or ComplexF64")
	end

	local bctype_code::Int = parse_bctype(bctype)

	if  isa(fs, Matrix{ComplexF64}) && bctype_code == 1
		error("Complex spline doesn't have natural spline. (bctype = 1/natural)")
	end

    spline = _spline_setup(xs, fs, Int32(bctype_code))
	
	finalizer(_destroy_spline, spline)
    
	return spline
end


function _spline_eval(spline::RealSplineType, x::Float64, derivs::Int=0)
	# x -> Float64
	# Returns a vector of Float64 (nqty)
	if derivs == 0
		f = Vector{Float64}(undef, spline.nqty)

		ccall((:spline_c_eval, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{Float64}),
			spline.handle, x, f)
		return f
	elseif derivs == 1
		f = Vector{Float64}(undef, spline.nqty)
		f1 = Vector{Float64}(undef, spline.nqty)
		ccall((:spline_c_eval_deriv, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}),
			spline.handle, x, f, f1)
		return f, f1
	elseif derivs == 2
		f = Vector{Float64}(undef, spline.nqty)
		f1 = Vector{Float64}(undef, spline.nqty)
		f2 = Vector{Float64}(undef, spline.nqty)
		ccall((:spline_c_eval_deriv2, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
			spline.handle, x, f, f1, f2)
		return f, f1, f2
	elseif derivs == 3
		f = Vector{Float64}(undef, spline.nqty)
		f1 = Vector{Float64}(undef, spline.nqty)
		f2 = Vector{Float64}(undef, spline.nqty)
		f3 = Vector{Float64}(undef, spline.nqty)
		ccall((:spline_c_eval_deriv3, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
			spline.handle, x, f, f1, f2, f3)
		return f, f1, f2, f3
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3.")
	end
end

function _spline_eval(spline::RealSplineType, xs::Vector{Float64}, derivs::Int=0)
    # xs -> Float64 (any length)
    # Returns a matrix of Float64 (length(xs), nqty)
    n = length(xs)
    fs = Matrix{Float64}(undef, n, spline.nqty)
	f = Vector{Float64}(undef, spline.nqty)
	if derivs > 0
		fs1 = Matrix{Float64}(undef, n, spline.nqty)
		f1 = Vector{Float64}(undef, spline.nqty)
	end
	if derivs > 1
		fs2 = Matrix{Float64}(undef, n, spline.nqty)
		f2 = Vector{Float64}(undef, spline.nqty)
	end
	if derivs > 2
		fs3 = Matrix{Float64}(undef, n, spline.nqty)
		f3 = Vector{Float64}(undef, spline.nqty)
	end
	for i in 1:n
		if derivs == 0
			f = _spline_eval(spline, xs[i], 0)
			fs[i, :] = f
		elseif derivs == 1
			f, f1 = _spline_eval(spline, xs[i], 1)
			fs[i, :] = f
			fs1[i, :] = f1
		elseif derivs == 2
			f, f1, f2 = _spline_eval(spline, xs[i], 2)
			fs[i, :] = f
			fs1[i, :] = f1
			fs2[i, :] = f2
		elseif derivs == 3
			f, f1, f2, f3 = _spline_eval(spline, xs[i], 3)
			fs[i, :] = f
			fs1[i, :] = f1
			fs2[i, :] = f2
			fs3[i, :] = f3
		else
			error("Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3.")
		end
	end
	if derivs == 0
		return fs
	elseif derivs == 1
		return fs, fs1
	elseif derivs == 2
		return fs, fs1, fs2
	elseif derivs == 3
		return fs, fs1, fs2, fs3
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3.")
	end
end


function _spline_eval(spline::ComplexSplineType, x::Float64, derivs::Int=0)
	# x -> Float64
	# Returns a vector of ComplexF64 (nqty)
	f1 = Vector{ComplexF64}(undef, spline.nqty)
	f2 = Vector{ComplexF64}(undef, spline.nqty)
	f3 = Vector{ComplexF64}(undef, spline.nqty)
	if derivs == 0
		f = Vector{ComplexF64}(undef, spline.nqty)
		ccall((:cspline_c_eval, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{ComplexF64}),
			spline.handle, x, f)
		return f
	elseif derivs == 1
		f = Vector{ComplexF64}(undef, spline.nqty)
		f1 = Vector{ComplexF64}(undef, spline.nqty)
		ccall((:cspline_c_eval_deriv, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{ComplexF64}, Ptr{ComplexF64}),
			spline.handle, x, f, f1)
		return f, f1
	elseif derivs == 2
		f = Vector{ComplexF64}(undef, spline.nqty)
		f1 = Vector{ComplexF64}(undef, spline.nqty)
		f2 = Vector{ComplexF64}(undef, spline.nqty)
		ccall((:cspline_c_eval_deriv2, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}),
			spline.handle, x, f, f1, f2)
		return f, f1, f2
	elseif derivs == 3
		f = Vector{ComplexF64}(undef, spline.nqty)
		f1 = Vector{ComplexF64}(undef, spline.nqty)
		f2 = Vector{ComplexF64}(undef, spline.nqty)
		f3 = Vector{ComplexF64}(undef, spline.nqty)
		ccall((:cspline_c_eval_deriv3, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}),
			spline.handle, x, f, f1, f2, f3)
		return f, f1, f2, f3
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3.")
	end
end

function _spline_eval(spline::ComplexSplineType, xs::Vector{Float64}, derivs::Int=0)
    # xs -> Float64 (any length)
    # Returns a matrix of ComplexF64 (length(xs), nqty)
    n = length(xs)
    fs = Matrix{ComplexF64}(undef, n, spline.nqty)
    f = Vector{ComplexF64}(undef, spline.nqty)
	if derivs > 0
		fs1 = Matrix{ComplexF64}(undef, n, spline.nqty)
		f1 = Vector{ComplexF64}(undef, spline.nqty)
	end
	if derivs > 1
		fs2 = Matrix{ComplexF64}(undef, n, spline.nqty)
		f2 = Vector{ComplexF64}(undef, spline.nqty)
	end
	if derivs > 2
		fs3 = Matrix{ComplexF64}(undef, n, spline.nqty)
		f3 = Vector{ComplexF64}(undef, spline.nqty)
	end
    for i in 1:n
		if derivs == 0
			f = _spline_eval(spline, xs[i], 0)
			fs[i, :] = f
		elseif derivs == 1
			f, f1 = _spline_eval(spline, xs[i], 1)
			fs[i, :] = f
			fs1[i, :] = f1
		elseif derivs == 2
			f, f1, f2 = _spline_eval(spline, xs[i], 2)
			fs[i, :] = f
			fs1[i, :] = f1
			fs2[i, :] = f2
		elseif derivs == 3
			f, f1, f2, f3 = _spline_eval(spline, xs[i], 3)
			fs[i, :] = f
			fs1[i, :] = f1
			fs2[i, :] = f2
			fs3[i, :] = f3
		else
			error("Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3.")
		end
	end
	if derivs == 0
		return fs
	elseif derivs == 1
		return fs, fs1
	elseif derivs == 2
		return fs, fs1, fs2
	elseif derivs == 3
		return fs, fs1, fs2, fs3
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3.")
	end
end




#internal version




function spline_eval(spline::CubicSplineType, x, derivs::Int=0)
	"""
	# spline_eval(spline, x)
	## Arguments:
	- `spline`: A `Spline` object created by `spline_setup`.
	- `x`: A Float64 value or a vector of Float64 values representing the x-coordinates to evaluate the spline at.
	## Returns:
	- If `x` is a single Float64 value, returns a vector of Float64 values representing the function values at that x-coordinate.
	- If `x` is a vector of Float64 values, returns a matrix of Float64 values where each row corresponds to the function values at
	  the respective x-coordinate in `x`.
	- Depending on the derivatives requested, it may return additional vectors for the first, second, or third derivatives.
	"""
	return _spline_eval(spline, x, derivs)
end



function _spline_integrate!(spline::RealSplineType)
    spline._fsi = Matrix{Float64}(undef, spline.mx + 1, spline.nqty)

    ccall((:spline_c_int, libspline), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}),
          spline.handle, spline._fsi)

    return
end

function _spline_integrate!(spline::ComplexSplineType)
    spline._fsi = Matrix{ComplexF64}(undef, spline.mx + 1, spline.nqty)


    ccall((:cspline_c_int, libspline), Cvoid,
          (Ptr{Cvoid}, Ptr{ComplexF64}), 
          spline.handle, spline._fsi)

    return
end

function spline_integrate!(spline::CubicSplineType)
    """
    spline_integrate!(spline)

	## Arguments:
	- `spline`: A mutable `Spline` object".

	## Returns:
	- Nothing. Updates `spline._fsi` in place so that  
	`spline._fsi[i, :]` equals `∫_{xs[1]}^{xs[i]} f(x) dx` for each component.
	"""

    _spline_integrate!(spline)
end

end