module BicubicSpline

const libdir = joinpath(@__DIR__, "..", "..", "deps")
const libspline = joinpath(libdir, "libspline")

export bicube_setup, bicube_eval

mutable struct BicubicSplineType
    handle::Ptr{Cvoid}
    xs::Vector{Float64}
    ys::Vector{Float64}
    fs::Array{Float64,3} # 3D array for bicubic spline values
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


function _bicube_setup(xs::Vector{Float64}, ys::Vector{Float64}, fs::Array{Float64,3}, bctypex::Int32, bctypey::Int32)
    # xs -> Float64 (mx)
    # ys -> Float64 (my)
    # fs -> Float64 (mx, my, nqty)
    if length(xs) != size(fs, 1) || length(ys) != size(fs, 2)
        error("Length of xs must match number of rows in fs and length of ys must match number of columns in fs")
    end
    mx = length(xs) - 1
    my = length(ys) - 1
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

"""
    bicube_setup(xs, ys, fs, bctypex=1, bctypey=1)

Set up a bicubic spline interpolation for 2D data.

# Arguments

  - `xs`: Vector of Float64 values representing the x-coordinates (must be monotonically increasing)
  - `ys`: Vector of Float64 values representing the y-coordinates (must be monotonically increasing)
  - `fs`: 3D array of Float64 values with dimensions (nx, ny, nqty) representing function values
  - `bctypex`: Boundary condition type for x-direction (default is 1)
  - `bctypey`: Boundary condition type for y-direction (default is 1)

      + 1: Natural spline (zero second derivative at boundaries)
      + 2: Periodic spline
      + 3: Extrapolated spline
      + 4: Not-a-knot spline

# Returns

  - A `BicubicSplineType` object that can be used for 2D evaluation

# Examples

```julia
# Create 2D grid
xs = collect(range(0.0; stop=2π, length=20))
ys = collect(range(0.0; stop=2π, length=20))

# Create 3D function data array
fs = zeros(20, 20, 1)
for i in 1:20, j in 1:20
    fs[i, j, 1] = sin(xs[i]) * cos(ys[j])
end

# Set up bicubic spline
bcspline = bicube_setup(xs, ys, fs, 1, 1)
```
"""
function bicube_setup(xs, ys, fs, bctypex::Int=1, bctypey::Int=1)
    if !isa(xs, Vector{Float64}) || !isa(ys, Vector{Float64}) || !isa(fs, Array{Float64,3})
        error("xs must be a vector of Float64, ys must be a vector of Float64, and fs must be a 3D array of Float64")
    end
    bicube = _bicube_setup(xs, ys, fs, Int32(bctypex), Int32(bctypey))
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

"""
    bicube_eval(bicube, x, y, derivs=0)

Evaluate a bicubic spline at given 2D points.

# Arguments

  - `bicube`: A `BicubicSplineType` object created by `bicube_setup`
  - `x`: Float64 value or vector of x-coordinates to evaluate
  - `y`: Float64 value or vector of y-coordinates to evaluate
  - `derivs`: Integer specifying derivative level (default is 0):

      + 0: Function values only
      + 1: Function values and first derivatives (∂f/∂x, ∂f/∂y)
      + 2: Function values, first and second derivatives (∂f/∂x, ∂f/∂y, ∂²f/∂x², ∂²f/∂x∂y, ∂²f/∂y²)

# Returns

  - If `derivs=0`: 3D array of function values (length(x) × length(y) × nqty)
  - If `derivs=1`: Tuple of (values, x_derivatives, y_derivatives)
  - If `derivs=2`: Tuple of (values, x_derivatives, y_derivatives, xx_derivatives, xy_derivatives, yy_derivatives)

# Examples

```julia
# Evaluate at single point
bcspline = bicube_setup(xs, ys, fs, 1, 1)
f_vals = bicube_eval(bcspline, π / 2, π / 4)

# Evaluate on grid
x_eval = collect(range(0, 2π; length=50))
y_eval = collect(range(0, 2π; length=50))
f_vals = bicube_eval(bcspline, x_eval, y_eval)

# Get derivatives
f_vals, fx, fy = bicube_eval(bcspline, x_eval, y_eval, 1)
```
"""
function bicube_eval(bicube::BicubicSplineType, x, y, derivs::Int=0)
    return _bicube_eval(bicube, x, y, derivs)
end

end # module BicubicSpline