# fortran function names
call_spline_c_create(::Type{Float64}, mx, nqty, h)    = ccall((:spline_c_create,  libspline), Cvoid, (Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)
call_spline_c_create(::Type{ComplexF64}, mx, nqty, h) = ccall((:cspline_c_create, libspline), Cvoid, (Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)
call_spline_c_setup(::Type{Float64}, spline, xs, fs)    = ccall((:spline_c_setup,  libspline), Cvoid, (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}), spline.handle, xs, fs)
call_spline_c_setup(::Type{ComplexF64}, spline, xs, fs) = ccall((:cspline_c_setup, libspline), Cvoid, (Ptr{Cvoid}, Ptr{Float64}, Ptr{ComplexF64}), spline.handle, xs, fs)
call_spline_c_fit(::Type{Float64}, spline)    = ccall((:spline_c_fit,  libspline), Cvoid, (Ptr{Cvoid}, Int32, Ptr{Float64}), spline.handle, spline.bctype, spline._fs1)
call_spline_c_fit(::Type{ComplexF64}, spline) = ccall((:cspline_c_fit, libspline), Cvoid, (Ptr{Cvoid}, Int32, Ptr{ComplexF64}), spline.handle, spline.bctype, spline._fs1)
call_spline_c_eval(::Type{Float64}, spline, x, f) = ccall((:spline_c_eval, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{Float64}), spline.handle, x, f)
call_spline_c_eval(::Type{ComplexF64}, spline, x, f) = ccall((:cspline_c_eval, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{ComplexF64}), spline.handle, x, f)
call_spline_c_eval(::Type{Float64}, spline, x, f, f1) = ccall((:spline_c_eval_deriv, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}), spline.handle, x, f, f1)
call_spline_c_eval(::Type{ComplexF64}, spline, x, f, f1) = ccall((:cspline_c_eval_deriv, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{ComplexF64}, Ptr{ComplexF64}), spline.handle, x, f, f1)
call_spline_c_eval(::Type{Float64}, spline, x, f, f1, f2) = ccall((:spline_c_eval_deriv2, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), spline.handle, x, f, f1, f2)
call_spline_c_eval(::Type{ComplexF64}, spline, x, f, f1, f2) = ccall((:cspline_c_eval_deriv2, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), spline.handle, x, f, f1, f2)
call_spline_c_eval(::Type{Float64}, spline, x, f, f1, f2, f3) = ccall((:spline_c_eval_deriv3, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), spline.handle, x, f, f1, f2, f3)
call_spline_c_eval(::Type{ComplexF64}, spline, x, f, f1, f2, f3) = ccall((:cspline_c_eval_deriv3, libspline), Cvoid, (Ptr{Cvoid}, Float64, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), spline.handle, x, f, f1, f2, f3)
call_spline_c_int(::Type{Float64}, spline) = ccall((:spline_c_int, libspline), Cvoid, (Ptr{Cvoid}, Ptr{Float64}), spline.handle, spline._fsi)
call_spline_c_int(::Type{ComplexF64}, spline) = ccall((:cspline_c_int, libspline), Cvoid, (Ptr{Cvoid}, Ptr{ComplexF64}), spline.handle, spline._fsi)
call_spline_c_destroy(::Type{Float64}, spline)    = ccall((:spline_c_destroy,  libspline), Cvoid, (Ptr{Cvoid},), spline.handle)
call_spline_c_destroy(::Type{ComplexF64}, spline) = ccall((:cspline_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), spline.handle)

mutable struct CubicSpline{T<:Union{Float64, ComplexF64}}
	handle::Ptr{Cvoid}
	_xs::Vector{Float64}
	_fs::Matrix{T}
	mx::Int
	nqty::Int
    bctype::Int32  # Boundary condition type
	_fsi::Matrix{T} # To store integrals at gridpoint
	_fs1::Matrix{T} # To store 1-deriv at gridpoint
end

function CubicSpline(unmanaged_handle::Ptr{Cvoid}, xs::Vector{Float64}, fs::Matrix{T}, mx::Int, nqty::Int) where {T<:Union{Float64, ComplexF64}}
	fsi = Matrix{T}(undef, 0, 0)
	fs1 = Matrix{T}(undef, 0, 0)
	return CubicSpline(unmanaged_handle, xs, fs, mx, nqty, zero(Int32), fsi, fs1)
end

@expose_fields CubicSpline xs fs fsi fs1

function _destroy_spline(spline::CubicSpline{T}) where {T<:Union{Float64, ComplexF64}}
    if spline.handle != C_NULL
        call_spline_c_destroy(T, spline)
		setfield!(spline, :handle, C_NULL)
    end
end

function CubicSpline(xs::Vector{Float64}, fs::Vector{<:Union{Float64, ComplexF64}}, bctype::Int32)
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

"""
spline_eval(spline::CubicSpline{T}, x, derivs::Int=0) where {T<:Union{Float64, ComplexF64}}
	## Arguments:
	- `spline`: A `Spline` object created by `spline_setup`.
	- `x`: A Float64 value or a vector of Float64 values representing the x-coordinates to evaluate the spline at.
	## Returns:
	- If `x` is a single Float64 value, returns a vector of Float64 values representing the function values at that x-coordinate.
	- If `x` is a vector of Float64 values, returns a matrix of Float64 values where each row corresponds to the function values at
	  the respective x-coordinate in `x`.
	- Depending on the derivatives requested, it may return additional vectors for the first, second, or third derivatives.
"""
function spline_eval(spline::CubicSpline{T}, x::Float64, derivs::Int=0) where {T<:Union{Float64, ComplexF64}}
	# x -> Float64
	# Returns a vector of T (nqty)
	@assert (derivs in 0:3) "Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3."

	f = Vector{T}(undef, spline.nqty)
	if derivs == 0
		call_spline_c_eval(T, spline, x, f)
		return f
	elseif derivs == 1
		f1 = similar(f)
		call_spline_c_eval(T, spline, x, f, f1)
		return f, f1
	elseif derivs == 2
		f1, f2 = similar(f), similar(f)
		call_spline_c_eval(T, spline, x, f, f1, f2)
		return f, f1, f2
	elseif derivs == 3
		f1, f2, f3 = similar(f), similar(f), similar(f)
		call_spline_c_eval(T, spline, x, f, f1, f2, f3)
		return f, f1, f2, f3
	end
end

function spline_eval(spline::CubicSpline{T}, xs::Vector{Float64}, derivs::Int=0) where {T<:Union{Float64, ComplexF64}}
    # xs -> Float64 (any length)
    # Returns a matrix of T (length(xs), nqty)
	@assert (derivs in 0:3) "Invalid number of derivatives requested: $derivs. Must be 0, 1, 2, or 3."

    n = length(xs)
    fs = Matrix{T}(undef, n, spline.nqty)
    f = Vector{T}(undef, spline.nqty)
	if derivs > 0
		fs1, f1 = similar(fs), similar(f)
	end
	if derivs > 1
		fs2, f2 = similar(fs), similar(f)
	end
	if derivs > 2
		fs3, f3 = similar(fs), similar(f)
	end
    for (i, x) in enumerate(xs)
		if derivs == 0
			call_spline_c_eval(T, spline, x, f)
			fs[i, :] .= f
		elseif derivs == 1
			call_spline_c_eval(T, spline, x, f, f1)
			fs[i, :] .= f
			fs1[i, :] .= f1
		elseif derivs == 2
			call_spline_c_eval(T, spline, x, f, f1, f2)
			fs[i, :] .= f
			fs1[i, :] .= f1
			fs2[i, :] .= f2
		elseif derivs == 3
			call_spline_c_eval(T, spline, x, f, f1, f2, f3)
			fs[i, :] .= f
			fs1[i, :] .= f1
			fs2[i, :] .= f2
			fs3[i, :] .= f3
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
	end
end

"""
    spline_integrate!(spline::CubicSpline{T}) where {T<:Union{Float64, ComplexF64}}

	## Arguments:
	- `spline`: A mutable `Spline` object".

	## Returns:
	- Nothing. Updates `spline._fsi` in place so that
	`spline._fsi[i, :]` equals `∫_{xs[1]}^{xs[i]} f(x) dx` for each component.
"""
spline_integrate!(spline::CubicSpline{T}) where {T<:Union{Float64, ComplexF64}} = call_spline_c_int(T, spline)
