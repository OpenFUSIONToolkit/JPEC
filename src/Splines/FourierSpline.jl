module FourierSpline

const libdir = joinpath(@__DIR__, "..", "..", "deps")
const libspline = joinpath(libdir, "libspline")

export fspline_setup, fspline_eval

mutable struct FourierSplineType
    handle::Ptr{Cvoid}
    xs::Vector{Float64}
    ys::Vector{Float64}
    fs::Array{Float64, 3}
    mx::Int64
    my::Int64
    mband::Int64
    nqty::Int64
    bctype::Int32
    fit_method::Int32

    # output save
    f::Vector{Float64}
    fx::Vector{Float64}
    fy::Vector{Float64}
    fxx::Vector{Float64}
    fxy::Vector{Float64}
    fyy::Vector{Float64}
end

function _destroy_fspline(fourier::FourierSplineType)
    if fourier.handle != C_NULL
        ccall((:fspline_c_destroy, libspline), Cvoid, (Ptr{Cvoid},), fourier.handle)
        fourier.handle = C_NULL
    end
end

function _MakeFourierSpline(mx::Int, my::Int, mband::Int, nqty::Int)
    h = Ref{Ptr{Cvoid}}()
    ccall((:fspline_c_create, libspline), Cvoid,
          (Int64, Int64, Int64, Int64, Ref{Ptr{Cvoid}}),
          mx, my, mband, nqty, h)

    handle = h[]
    if handle == C_NULL
        error("Failed to create fspline handle in Fortran library.")
    end


    f = Vector{Float64}(undef, nqty)
    fx = Vector{Float64}(undef, nqty)
    fy = Vector{Float64}(undef, nqty)
    fxx = Vector{Float64}(undef, nqty)
    fxy = Vector{Float64}(undef, nqty)
    fyy = Vector{Float64}(undef, nqty)

    # Return a partially initialized struct with empty data arrays
    return FourierSplineType(handle,
                             Vector{Float64}(undef, 0),
                             Vector{Float64}(undef, 0),
                             Array{Float64, 3}(undef, 0, 0, 0),
                             mx, my, mband, nqty, 0, 0,
                             f, fx, fy, fxx, fxy, fyy)
end

function _fspline_setup(xs::Vector{Float64}, ys::Vector{Float64}, fs::Array{Float64, 3}
    , mband::Int, bctype::Int32, fit_method::Int32, fit_flag::Bool)
    
    
    mx = length(xs) - 1
    my = length(ys) - 1
    nqty = size(fs, 3)
    fourier = _MakeFourierSpline(mx, my, mband, nqty)
    fourier.xs = xs
    fourier.ys = ys
    fourier.fs = fs
    fourier.bctype = Int32(bctype)
    fourier.fit_method = Int32(fit_method)

    ccall((:fspline_c_setup, libspline), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          fourier.handle, xs, ys, fs)
    if fit_method == 1
        ccall((:fspline_c_fit_1, libspline), Cvoid,
              (Ptr{Cvoid}, Int32, Bool),
              fourier.handle, bctype, fit_flag)
    elseif fit_method == 2
        ccall((:fspline_c_fit_2, libspline), Cvoid,
              (Ptr{Cvoid}, Int32, Bool),
              fourier.handle, bctype, fit_flag)
    else
        error("Internal error: Invalid fit_method passed to _fspline_setup.")
    end

    return fourier
end

function fspline_setup(xs::Vector{Float64}, ys::Vector{Float64}, fs::Array{Float64, 3}
    , mband::Int; bctype::Int=4, fit_method::Int=1, fit_flag::Bool=true)
    """
    # fspline_setup(xs, ys, fs, mband; bctype=1, fit_method=1)

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

    if bctype == 1
        error("Fourier spline  doesn't have natural spline. (bctype = 1)")
    end
    # Call the internal setup function that creates and fits the spline
    fourier = _fspline_setup(xs, ys, fs, mband, Int32(bctype), Int32(fit_method),fit_flag)

    # Add a finalizer to ensure the Fortran object is deallocated when the Julia object is garbage collected.
    finalizer(_destroy_fspline, fourier)

    return fourier
end

function _fspline_eval(fourier::FourierSplineType, x::Float64, y::Float64, derivs::Int)
    if derivs == 0
        f = Vector{Float64}(undef, fourier.nqty)
        ccall((:fspline_c_eval, libspline), Cvoid,
              (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}),
              fourier.handle, x, y, f)
        return f
    elseif derivs == 1
        f = Vector{Float64}(undef, fourier.nqty)
        fx = Vector{Float64}(undef, fourier.nqty)
        fy = Vector{Float64}(undef, fourier.nqty)
        ccall((:fspline_c_eval_deriv, libspline), Cvoid,
              (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
              fourier.handle, x, y, f, fx, fy)
        return f, fx, fy
    elseif derivs == 2
        f = Vector{Float64}(undef, fourier.nqty)
        fx = Vector{Float64}(undef, fourier.nqty)
        fy = Vector{Float64}(undef, fourier.nqty)
        fxx = Vector{Float64}(undef, fourier.nqty)
        fxy = Vector{Float64}(undef, fourier.nqty)
        fyy = Vector{Float64}(undef, fourier.nqty)
        ccall((:fspline_c_eval_deriv2, libspline), Cvoid,
              (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
              fourier.handle, x, y, f, fx, fy, fxx, fxy, fyy)
        return f, fx, fy, fxx, fxy, fyy
    else
        error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
    end
end






function _fspline_eval(fourier::FourierSplineType, xs::Vector{Float64}, ys::Vector{Float64}, derivs::Int)
    nx = length(xs)
    ny = length(ys)

    # Pre-allocate output arrays
    f_grid = Array{Float64, 3}(undef, nx, ny, fourier.nqty)
    if derivs > 0
        fx_grid = Array{Float64, 3}(undef, nx, ny, fourier.nqty)
        fy_grid = Array{Float64, 3}(undef, nx, ny, fourier.nqty)
    end
    if derivs > 1
        fxx_grid = Array{Float64, 3}(undef, nx, ny, fourier.nqty)
        fxy_grid = Array{Float64, 3}(undef, nx, ny, fourier.nqty)
        fyy_grid = Array{Float64, 3}(undef, nx, ny, fourier.nqty)
    end

    for j in 1:ny, i in 1:nx
        if derivs == 0
            f = _fspline_eval(fourier, xs[i], ys[j], 0)
            f_grid[i, j, :] = f
        elseif derivs == 1
            f, fx, fy = _fspline_eval(fourier, xs[i], ys[j], 1)
            f_grid[i, j, :] = f
            fx_grid[i, j, :] = fx
            fy_grid[i, j, :] = fy
        elseif derivs == 2
            f, fx, fy, fxx, fxy, fyy = _fspline_eval(fourier, xs[i], ys[j], 2)
            f_grid[i, j, :] = f
            fx_grid[i, j, :] = fx
            fy_grid[i, j, :] = fy
            fxx_grid[i, j, :] = fxx
            fxy_grid[i, j, :] = fxy
            fyy_grid[i, j, :] = fyy
        end
    end

    if derivs == 0
        return f_grid
    elseif derivs == 1
        return f_grid, fx_grid, fy_grid
    elseif derivs == 2
        return f_grid, fx_grid, fy_grid, fxx_grid, fxy_grid, fyy_grid
    else
        # This branch is technically unreachable due to checks in the public function
        error("Internal error during evaluation.")
    end
end

function _fspline_eval!(fourier::FourierSplineType, x::Float64, y::Float64, derivs::Int=0)
    # Modifies fourier.f, .fx, .fy, etc. in place.
    if derivs == 0
        ccall((:fspline_c_eval, libspline), Cvoid,
              (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}),
              fourier.handle, x, y, fourier.f)
    elseif derivs == 1
        ccall((:fspline_c_eval_deriv, libspline), Cvoid,
              (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
              fourier.handle, x, y, fourier.f, fourier.fx, fourier.fy)
    elseif derivs == 2
        ccall((:fspline_c_eval_deriv2, libspline), Cvoid,
              (Ptr{Cvoid}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
              fourier.handle, x, y, fourier.f, fourier.fx, fourier.fy, fourier.fxx, fourier.fxy, fourier.fyy)
    else
        error("Invalid number of derivatives requested: $derivs. Must be 0, 1, or 2.")
    end
    return
end

function fspline_eval(fourier::FourierSplineType, x, y; derivs::Int=0)
    """
    # fspline_eval(spl, x, y; derivs=0)

    Evaluates a fitted Fourier-Spline at given coordinates.

    ## Arguments:
    - `spl`: A `FourierSplineType` object from `fspline_setup`.
    - `x`: A `Float64` or a `Vector{Float64}` of x-coordinates.
    - `y`: A `Float64` or a `Vector{Float64}` of y-coordinates.

    ## Keyword Arguments:
    - `derivs`: Number of derivatives to compute.
        - 0: Function value `f` (default).
        - 1: `f`, `fx`, `fy`.
        - 2: `f`, `fx`, `fy`, `fxx`, `fxy`, `fyy`.

    ## Returns:
    - If `x`, `y` are scalars: A tuple containing the function value(s) and any requested derivatives. If nqty=1, the results are scalars, otherwise they are vectors.
    - If `x`, `y` are vectors: A tuple of 3D arrays for the function values and derivatives on the grid defined by `x` and `y`.
    """
    if derivs < 0 || derivs > 2
        error("Keyword `derivs` must be 0, 1, or 2.")
    end

    # Handle evaluation for single points or grids
    results = _fspline_eval(fourier, x, y, derivs)
    
    # Squeeze the result if nqty is 1 and inputs were scalars
    if fourier.nqty == 1 && isa(x, Real) && isa(y, Real)
        if isa(results, Tuple)
            return map(vec -> vec[1], results)
        else
            return results[1]
        end
    end

    return results
end


function fspline_eval!(fourier::FourierSplineType, x::Float64, y::Float64, derivs::Int=0)
    """
    # fspline_eval!(fourier, x, y; derivs=0)

    Evaluates the Fourier-spline at a point (x, y) and stores the result
    in-place in the `fourier` object's fields (`.f`, `.fx`, etc.).
    This method avoids memory allocation.

    ## Arguments:
    - `fourier`: A `FourierSplineType` object from `fspline_setup`.
    - `x`: A `Float64` x-coordinate.
    - `y`: A `Float64` y-coordinate.
    - `derivs`: Number of derivatives to compute (0, 1, or 2).
    """
    if derivs < 0 || derivs > 2
        error("Keyword `derivs` must be 0, 1, or 2.")
    end
    _fspline_eval!(fourier, x, y, derivs)
end



end # module FourierSpline