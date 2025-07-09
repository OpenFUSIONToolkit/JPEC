module BicubicSpline

const libdir = joinpath(@__DIR__, "..", "..", "deps")
const libspline = joinpath(libdir, "libspline")

export bicube_setup, bicube_eval

mutable struct BicubicSplineType
	handle::Ptr{Cvoid}
	xs::Vector{Float64}
	ys::Vector{Float64}
	fs::Array{Float64, 3} # 3D array for bicubic spline values
	mx::Int64
	my::Int64
	nqty::Int64
	ix::Int32  # Index of x position in the spline
	iy::Int32  # Index of y position in the spline
	bctypex::Int32  # Boundary condition type for x
	bctypey::Int32  # Boundary condition type for y


	fsx::Array{Float64, 3}
	fsy::Array{Float64, 3}
	fsxy::Array{Float64, 3}


	#output save
	f::Vector{Float64}
    fx::Vector{Float64}
    fy::Vector{Float64}
    fxx::Vector{Float64}
	fxy::Vector{Float64}
	fyy::Vector{Float64}
end

function _destroy_bicubic_spline(bicube::BicubicSplineType)
    if bicube.handle != C_NULL
        ccall((:bicube_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), bicube.handle)
        bicube.handle = C_NULL
    end
end

function _MakeBicubicSpline(mx::Int64, my::Int64, nqty::Int64)
	h = Ref{Ptr{Cvoid}}()
	ccall((:bicube_c_create, libspline), Cvoid,
		(Int64, Int64, Int64, Ref{Ptr{Cvoid}}), mx, my, nqty, h)

	fsx = Array{Float64, 3}(undef, 0,0,0)
	fsy = Array{Float64, 3}(undef, 0,0,0)
	fsxy = Array{Float64, 3}(undef, 0,0,0)

	f = Vector{Float64}(undef, nqty)
	fx = Vector{Float64}(undef, nqty)
	fy = Vector{Float64}(undef, nqty)
	fxx = Vector{Float64}(undef, nqty)
	fxy = Vector{Float64}(undef, nqty)
	fyy = Vector{Float64}(undef, nqty)


	return BicubicSplineType(h[], Vector{Float64}(undef, 0), Vector{Float64}(undef, 0), 
		Array{Float64,3}(undef, 0, 0, 0), mx, my, nqty, 0, 0, 0, 0
		,fsx,fsy, fsxy, f, fx, fy, fxx, fxy, fyy)
end




function _bicube_setup(xs::Vector{Float64}, ys::Vector{Float64}, fs::Array{Float64, 3}, bctypex::Int32, bctypey::Int32)
	# xs -> Float64 (mx)
	# ys -> Float64 (my)
	# fs -> Float64 (mx, my, nqty)
	if length(xs) != size(fs, 1) || length(ys) != size(fs, 2)
		error("Length of xs must match number of rows in fs and length of ys must match number of columns in fs")
	end
	mx = length(xs)-1
	my = length(ys)-1
	nqty = size(fs, 3)
	bicube = _MakeBicubicSpline(mx, my, nqty)
	bicube.xs = xs
	bicube.ys = ys
	bicube.fs = fs
	bicube.bctypex = Int32(bctypex)
	bicube.bctypey = Int32(bctypey)

	bicube.fsx = Array{Float64, 3}(undef, mx+1, my+1, nqty)
	bicube.fsy = Array{Float64, 3}(undef, mx+1, my+1, nqty)
	bicube.fsxy = Array{Float64, 3}(undef, mx+1, my+1, nqty)


	ccall((:bicube_c_setup, libspline), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		bicube.handle, xs, ys, fs)

	ccall((:bicube_c_fit, libspline), Cvoid,
		(Ptr{Cvoid}, Int32, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
		bicube.handle, bctypex, bctypey, bicube.fsx, bicube.fsy, bicube.fsxy)
	
	
	
	
	return bicube
end

function bicube_setup(xs, ys, fs, bctypex::Int=4, bctypey::Int=4)
	"""
	# bicube_setup(xs, ys, fs, bctypex=0, bctypey=0)
	## Arguments:
	- `xs`: A vector of Float64 values representing the x-coordinates.
	- `ys`: A vector of Float64 values representing the y-coordinates.
	- `fs`: A 3D array of Float64 values representing the function values at the (x,y) coordinates.
	- `bctypex`: An integer specifying the boundary condition type for x (Defaule 4, not a knot)
	- `bctypey`: An integer specifying the boundary condition type for y  (Defaule 4, not a knot) 
	## Returns:
	- A `BicubicSpline` object containing the spline handle, x-coordinates, y-coordinates,
	function values, number of x-coordinates, number of y-coordinates, number of quantities,
	and boundary condition types.
	"""
	if !isa(xs, Vector{Float64}) || !isa(ys, Vector{Float64}) || !isa(fs, Array{Float64, 3})
		error("xs must be a vector of Float64, ys must be a vector of Float64, and fs must be a 3D array of Float64")
	end
	bicube = _bicube_setup(xs, ys, fs, Int32(bctypex), Int32(bctypey))
	
	finalizer(_destroy_bicubic_spline, bicube)
	
	return bicube
end

function _bicube_eval(bicube::BicubicSplineType, x::Float64, y::Float64, derivs::Int=0)
	# x -> Float64
	# y -> Float64
	# Returns a vector of Float64 (nqty)
	f = Vector{Float64}(undef, bicube.nqty)
	fx = Vector{Float64}(undef, bicube.nqty)
	fy = Vector{Float64}(undef, bicube.nqty)
	fxx = Vector{Float64}(undef, bicube.nqty)
	fxy = Vector{Float64}(undef, bicube.nqty)
	fyy = Vector{Float64}(undef, bicube.nqty)
	if derivs == 0
		ccall((:bicube_c_eval, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Int32, Int32),
			bicube.handle, x, y, f, bicube.ix, bicube.iy)
	elseif derivs == 1
		ccall((:bicube_c_eval_deriv, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32, Int32),
			bicube.handle, x, y, f, fx, fy, bicube.ix, bicube.iy)
	elseif derivs == 2
		ccall((:bicube_c_eval_deriv2, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32, Int32),
			bicube.handle, x, y, f, fx, fy, fxx, fxy, fyy, bicube.ix, bicube.iy)
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
	end
	if derivs == 0
		return f
	elseif derivs == 1
		return f, fx, fy
	elseif derivs == 2
		return f, fx, fy, fxx, fxy, fyy
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
	end
end


function _bicube_eval(bicube::BicubicSplineType, xs::Vector{Float64}, ys::Vector{Float64}, derivs::Int=0)
	# xs -> Float64 (any length)
	# ys -> Float64 (any length)
	# Returns a matrix of Float64 (length(xs), length(ys), nqty)

	n = length(xs)
	m = length(ys)
	# fs = Array{Float64}(undef, n, m, bicube.nqty)
	# fsx = Array{Float64}(undef, n, m, bicube.nqty)
	# fsy = Array{Float64}(undef, n, m, bicube.nqty)
	# fsxx = Array{Float64}(undef, n, m, bicube.nqty)
	# fsxy = Array{Float64}(undef, n, m, bicube.nqty)
	# fsyy = Array{Float64}(undef, n, m, bicube.nqty)
	# f = Vector{Float64}(undef, bicube.nqty)
	# fx = Vector{Float64}(undef, bicube.nqty)
	# fy = Vector{Float64}(undef, bicube.nqty)
	# fxx = Vector{Float64}(undef, bicube.nqty)
	# fxy = Vector{Float64}(undef, bicube.nqty)
	# fyy = Vector{Float64}(undef, bicube.nqty)
	fs = Array{Float64}(undef, n, m, bicube.nqty)
	f = Vector{Float64}(undef, bicube.nqty)
	if derivs > 0
		fsx = Array{Float64}(undef, n, m, bicube.nqty)
		fsy = Array{Float64}(undef, n, m, bicube.nqty)
		fx = Vector{Float64}(undef, bicube.nqty)
		fy = Vector{Float64}(undef, bicube.nqty)
	end
	if derivs > 1
		fsxx = Array{Float64}(undef, n, m, bicube.nqty)
		fsxy = Array{Float64}(undef, n, m, bicube.nqty)
		fsyy = Array{Float64}(undef, n, m, bicube.nqty)
		fxx = Vector{Float64}(undef, bicube.nqty)
		fxy = Vector{Float64}(undef, bicube.nqty)
		fyy = Vector{Float64}(undef, bicube.nqty)
	end
	for i in 1:n
		for j in 1:m
			if derivs == 0
				f = _bicube_eval(bicube, xs[i], ys[j], 0)
				fs[i, j, :] = f
			elseif derivs == 1
				f, fx, fy = _bicube_eval(bicube, xs[i], ys[j], 1)
				fs[i, j, :] = f
				fsx[i, j, :] = fx
				fsy[i, j, :] = fy
			elseif derivs == 2
				f, fx, fy, fxx, fxy, fyy = _bicube_eval(bicube, xs[i], ys[j], 2)
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



function _bicube_eval!(bicube::BicubicSplineType, x::Float64, y::Float64, derivs::Int=0)
	# x -> Float64
	# y -> Float64
    # Modifies bicube.f, .fx, .fy, etc. in place.
	if derivs == 0
		ccall((:bicube_c_eval, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Int32, Int32),
			bicube.handle, x, y, bicube.f, bicube.ix, bicube.iy)
	elseif derivs == 1
		ccall((:bicube_c_eval_deriv, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32, Int32),
			bicube.handle, x, y, bicube.f, bicube.fx, bicube.fy, bicube.ix, bicube.iy)
	elseif derivs == 2
		ccall((:bicube_c_eval_deriv2, libspline), Cvoid,
			(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32, Int32),
			bicube.handle, x, y, bicube.f, bicube.fx, bicube.fy, bicube.fxx, bicube.fxy, bicube.fyy, bicube.ix, bicube.iy)
	else
		error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
	end
    return
end

function bicube_eval(bicube::BicubicSplineType, x, y, derivs::Int=0)
	"""
	# bicube_eval(bicube, x, y)
	## Arguments:
	- `bicube`: A `BicubicSpline` object created by `bicube_setup`.
	- `x`: A Float64 value or a vector of Float64 values representing the x-coordinates to evaluate the bicubic spline at.
	- `y`: A Float64 value or a vector of Float64 values representing the y-coordinates to evaluate the bicubic spline at.
	## Returns:
	- If `x` and `y` are single Float64 values, returns a vector of Float64 values representing the function values at that (x,y) coordinate.
	- If `x` and `y` are vectors of Float64 values, returns a 3D array of Float64 values where each slice corresponds to the function values at
	the respective (x,y) coordinates in `x` and `y`.
	"""
	return _bicube_eval(bicube, x, y, derivs)
end


function bicube_eval!(bicube::BicubicSplineType, x::Float64, y::Float64, derivs::Int=0)
    """
    # bicube_eval!(bicube, x, y; derivs=0)

    Evaluates the bicubic spline at a point (x, y) and stores the result
    in-place in the `bicube` object's fields (`.f`, `.fx`, etc.).
    This method avoids memory allocation and is suitable for performance-critical code.

    ## Arguments:
    - `bicube`: A `BicubicSpline` object from `bicube_setup`.
    - `x`: A `Float64` x-coordinate.
    - `y`: A `Float64` y-coordinate.
    - `derivs`: Number of derivatives to compute (0, 1, or 2).
    """
    _bicube_eval!(bicube, x, y, derivs)
end

end # module BicubicSpline