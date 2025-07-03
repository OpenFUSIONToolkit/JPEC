module SplinesMod

include("CubicSpline.jl")
include("BicubicSpline.jl")

using .CubicSpline: spline_setup, spline_eval
using .BicubicSpline: bicube_setup, bicube_eval

export spline_setup, spline_eval
export bicube_setup, bicube_eval

end