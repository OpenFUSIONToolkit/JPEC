module SplinesMod

const libdir = joinpath(@__DIR__, "..", "..", "deps")
const libspline = joinpath(libdir, "libspline")

include("Helper.jl")

include("CubicSpline.jl")
include("BicubicSpline.jl")
include("FourierSpline.jl")

export spline_setup, spline_eval, spline_integrate!, CubicSpline
export bicube_setup, bicube_eval, BicubicSpline
export fspline_setup, fspline_eval, FourierSpline

end
