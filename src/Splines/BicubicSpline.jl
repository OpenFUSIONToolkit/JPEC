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
end

function _MakeBicubicSpline(mx::Int64, my::Int64, nqty::Int64)
	h = Ref{Ptr{Cvoid}}()
	ccall((:bicube_c_create, libspline), Cvoid,
		(Int64, Int64, Int64, Ref{Ptr{Cvoid}}), mx, my, nqty, h)
	return BicubicSplineType(h[], Vector{Float64}(undef, mx), Vector{Float64}(undef, my), 
		Array{Float64,3}(undef, mx, my, nqty), mx, my, nqty, 0, 0, 0, 0)
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

	ccall((:bicube_c_setup, libspline), Cvoid,
		(Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		bicube.handle, xs, ys, fs)

	ccall((:bicube_c_fit, libspline), Cvoid,
		(Ptr{Cvoid}, Int32, Int32), bicube.handle, bctypex, bctypey)
	return bicube
end

function bicube_setup(xs, ys, fs, bctypex::Int=1, bctypey::Int=1)
	"""
	# bicube_setup(xs, ys, fs, bctypex=0, bctypey=0)
	## Arguments:
	- `xs`: A vector of Float64 values representing the x-coordinates.
	- `ys`: A vector of Float64 values representing the y-coordinates.
	- `fs`: A 3D array of Float64 values representing the function values at the (x,y) coordinates.
	- `bctypex`: An integer specifying the boundary condition type for x (default is 0).
	- `bctypey`: An integer specifying the boundary condition type for y (default is 0).
	## Returns:
	- A `BicubicSpline` object containing the spline handle, x-coordinates, y-coordinates,
	function values, number of x-coordinates, number of y-coordinates, number of quantities,
	and boundary condition types.
	"""
	if !isa(xs, Vector{Float64}) || !isa(ys, Vector{Float64}) || !isa(fs, Array{Float64, 3})
		error("xs must be a vector of Float64, ys must be a vector of Float64, and fs must be a 3D array of Float64")
	end
	bicube = _bicube_setup(xs, ys, fs, Int32(bctypex), Int32(bctypey))
	return bicube
end

function _bicube_eval(bicube::BicubicSplineType, x::Float64, y::Float64)
	# x -> Float64
	# y -> Float64
	# Returns a vector of Float64 (nqty)
	f = Vector{Float64}(undef, bicube.nqty)
	ccall((:bicube_c_eval, libspline), Cvoid,
		(Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Int32, Int32),
		bicube.handle, x, y, f, bicube.ix, bicube.iy)
	return f
end

function _bicube_eval(bicube::BicubicSplineType, xs::Vector{Float64}, ys::Vector{Float64})
	# xs -> Float64 (any length)
	# ys -> Float64 (any length)
	# Returns a matrix of Float64 (length(xs), length(ys), nqty)

	n = length(xs)
	m = length(ys)
	fs = Array{Float64}(undef, n, m, bicube.nqty)
	f = Vector{Float64}(undef, bicube.nqty)
	for i in 1:n
		for j in 1:m
			f = _bicube_eval(bicube, xs[i], ys[j])
			fs[i, j, :] = f
		end
	end
	return fs
end

function bicube_eval(bicube::BicubicSplineType, x, y)
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
	return _bicube_eval(bicube, x, y)
end

end # module BicubicSpline