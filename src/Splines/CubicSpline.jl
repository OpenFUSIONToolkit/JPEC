# fortran function names
call_spline_c_create(::Type{Float64}, mx, nqty, h)    = ccall((:spline_c_create,  libspline), Cvoid, (Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)
call_spline_c_create(::Type{ComplexF64}, mx, nqty, h) = ccall((:cspline_c_create, libspline), Cvoid, (Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)
call_spline_c_setup(::Type{Float64}, spline, xs, fs)    = ccall((:spline_c_setup,  libspline), Cvoid, (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}), spline.handle, xs, fs)
call_spline_c_setup(::Type{ComplexF64}, spline, xs, fs) = ccall((:cspline_c_setup, libspline), Cvoid, (Ptr{Cvoid}, Ptr{Float64}, Ptr{ComplexF64}), spline.handle, xs, fs)
call_spline_c_fit(::Type{Float64}, spline)    = ccall((:spline_c_fit,  libspline), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Float64}), spline.handle, spline.bctype, spline._fs1)
call_spline_c_fit(::Type{ComplexF64}, spline) = ccall((:cspline_c_fit, libspline), Cvoid, (Ptr{Cvoid}, Int32, Ptr{ComplexF64}), spline.handle, spline.bctype, spline._fs1)
call_spline_c_destroy(::Type{Float64}, spline)    = ccall((:spline_c_destroy,  libspline), Cvoid, (Ptr{Cvoid},), spline.handle)
call_spline_c_destroy(::Type{ComplexF64}, spline) = ccall((:cspline_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), spline.handle)

mutable struct CubicSpline{T<:Number}
	handle::Ptr{Cvoid}
	_xs::Vector{Float64}
	_fs::Matrix{T}
	mx::Int64
	nqty::Int64
    bctype::Int32  # Boundary condition type

	_fsi::Matrix{T} # To store integrals at gridpoint
	_fs1::Matrix{T} # To store 1-deriv at gridpoint
end

@expose_fields CubicSpline xs fs fsi fs1

function _destroy_spline(spline::CubicSpline{T}) where {T<:Number}
    if spline.handle != C_NULL
        call_spline_c_destroy(T, spline)
		setfield!(spline, :handle, C_NULL)
    end
end

function CubicSpline(xs::Vector{Float64}, fs::Vector{<:Number}, bctype::Int32)
	@assert length(xs) == length(fs) "Length of xs must match length of fs"
	fs_matrix = reshape(fs, length(xs), 1) # Convert to a column vector
	return CubicSpline(xs, fs_matrix, bctype)
end

function CubicSpline(xs::Vector{Float64}, fs::Matrix{T}, bctype::Int32) where {T <: Union{Float64, ComplexF64}}
	# xs -> Float64 (mx)
	# fs -> ComplexF64 (mx, nqty)
	@assert length(xs) == size(fs, 1) "Length of xs must match number of rows in fs"

	mx = length(xs)-1
	nqty = size(fs, 2)
	h = Ref{Ptr{Cvoid}}()
	call_spline_c_create(T, mx, nqty, h)

	fsi = Matrix{T}(undef, mx+1, nqty)
	fs1 = Matrix{T}(undef, mx+1, nqty)

	spline = CubicSpline(h[], xs, fs, mx, nqty, bctype, fsi, fs1)

	call_spline_c_setup(T, spline, xs, fs)
    call_spline_c_fit(T, spline)
	return spline
end

"""
    CubicSpline(xs, fs; bctype="not-a-knot")

    Arguments:
      - `xs`: A vector of Float64 values representing the x-coordinates.
      - `fs`: A vector or matrix of Float64/ComplexF64 values representing the function values at the x-coordinates.
    Keyword Arguments:
      - `bctype`: Boundary condition type for the cubic spline in `x`.
        - 1: Natural spline (default)
        - 2: Periodic spline
        - 3: Extrapolated spline
        - 4: "Not-a-knot" spline

	Returns:
      - A `CubicSpline` object containing the spline handle, x-coordinates, function values,
        number of x-coordinates, number of quantities, and index of x position in the spline.
"""
function CubicSpline(xs::Vector{Float64}, fs::Union{Vector{T}, Matrix{T}};
					 bctype::Union{String, Int}="not-a-knot") where {T <: Union{Float64, ComplexF64}}

	bctype_code = parse_bctype(bctype)

	if  T === ComplexF64 && bctype_code == 1
		error("Complex spline doesn't have natural spline. (bctype = 1/natural)")
	end

    spline = CubicSpline(xs, fs, Int32(bctype_code))

	finalizer(_destroy_spline, spline)

	return spline
end

# BCL 8/1/2025: UPDATED THROUGH HERE; THE BELOW WORKS BUT NEEDS TO BE JULIAFIED

function _spline_eval(spline::CubicSpline{Float64}, x::Float64, derivs::Int=0)
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

function _spline_eval(spline::CubicSpline{Float64}, xs::Vector{Float64}, derivs::Int=0)
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


function _spline_eval(spline::CubicSpline{ComplexF64}, x::Float64, derivs::Int=0)
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

function _spline_eval(spline::CubicSpline{ComplexF64}, xs::Vector{Float64}, derivs::Int=0)
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




function spline_eval(spline::CubicSpline, x, derivs::Int=0)
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



function _spline_integrate!(spline::CubicSpline{Float64})

    ccall((:spline_c_int, libspline), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}),
          spline.handle, spline._fsi)

    return
end

function _spline_integrate!(spline::CubicSpline{ComplexF64})

    ccall((:cspline_c_int, libspline), Cvoid,
          (Ptr{Cvoid}, Ptr{ComplexF64}),
          spline.handle, spline._fsi)

    return
end

function spline_integrate!(spline::CubicSpline)
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
