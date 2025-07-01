module GPECsplines

mutable struct Spline
	handle::Ptr{Cvoid}
	xs::Vector{Float64}
	fs::Matrix{Float64}
	mx::Int64
	nqty::Int64
	ix::Int32  # Index of x position in the spline
    bctype::Int32  # Boundary condition type
end

function _MakeSpline(mx::Int64, nqty::Int64)
	h = Ref{Ptr{Cvoid}}()
	ccall((:spline_c_create, "libspline"), Cvoid,
		(Int64, Int64, Ref{Ptr{Cvoid}}), mx, nqty, h)
	return Spline(h[], Vector{Float64}(undef, mx), Matrix{Float64}(undef, mx, nqty), mx, nqty, 0, 0)
end

function _spline_setup(xs::Vector{Float64}, fs::Vector{Float64}, bctype::Int32)
	# xs -> Float64 (mx)
	# fs -> Float64 (mx, nqty)
	if length(xs) != length(fs)
		error("Length of xs must match length of fs")
	end
	mx = length(xs)-1
	nqty = 1  # Default to 1 quantity if not specified
	spline = _MakeSpline(mx, nqty)
	spline.xs = xs
	# Convert fs to a matrix with one column
	fs_matrix = reshape(fs, mx, nqty)
	spline.fs = fs_matrix
    spline.bctype = bctype

	ccall((:spline_c_setup, "libspline"), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}),
		spline.handle, xs, fs_matrix)

    ccall((:spline_c_fit, "libspline"), Cvoid, 
        (Ptr{Cvoid}, Int32), spline.handle, bctype)
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
	spline = _MakeSpline(mx, nqty)
	spline.xs = xs
	spline.fs = fs

	ccall((:spline_c_setup, "libspline"), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}),
		spline.handle, xs, fs)

    ccall((:spline_c_fit, "libspline"), Cvoid, 
        (Ptr{Cvoid}, Int32), spline.handle, bctype)
	return spline
end

function spline_setup(xs, fs, bctype::Int=1)
    """
    # spline_setup(xs, fs, bctype=0)
    ## Arguments:
    - `xs`: A vector of Float64 values representing the x-coordinates.
    - `fs`: A vector or matrix of Float64 values representing the function values at the x-coordinates.
    - `bctype`: An integer specifying the boundary condition type (default is 0):
        - 1: Natural spline (default)
        - 2: Periodic spline
        - 3: Extrapolated spline
        - 4: not-a-knot spline
    ## Returns:
    - A `Spline` object containing the spline handle, x-coordinates, function values,
      number of x-coordinates, number of quantities, and index of x position in the spline.
    """
    spline = _spline_setup(xs, fs, Int32(bctype))
    return spline
end

function _spline_eval(spline::Spline, x::Float64)
	# x -> Float64
	# Returns a vector of Float64 (nqty)
	f = Vector{Float64}(undef, spline.nqty)
	ccall((:spline_c_eval, "libspline"), Cvoid,
		(Ptr{Cvoid}, Float64, Ptr{Float64}, Int32),
		spline.handle, x, f, spline.ix)
	return f
end

function _spline_eval(spline::Spline, xs::Vector{Float64})
    # xs -> Float64 (any length)
    # Returns a matrix of Float64 (length(xs), nqty)

    n = length(xs)
    fs = Matrix{Float64}(undef, n, spline.nqty)
    f = Vector{Float64}(undef, spline.nqty)
    for i in 1:n
        f = _spline_eval(spline, xs[i])
        fs[i, :] = f
    end
    return fs
end

function spline_eval(spline::Spline, x)
	"""
	# spline_eval(spline, x)
	## Arguments:
	- `spline`: A `Spline` object created by `spline_setup`.
	- `x`: A Float64 value or a vector of Float64 values representing the x-coordinates to evaluate the spline at.
	## Returns:
	- If `x` is a single Float64 value, returns a vector of Float64 values representing the function values at that x-coordinate.
	- If `x` is a vector of Float64 values, returns a matrix of Float64 values where each row corresponds to the function values at
	  the respective x-coordinate in `x`.
	"""
	return _spline_eval(spline, x)
end

end
