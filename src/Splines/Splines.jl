module SplinesMod

include("CubicSpline.jl")
include("BicubicSpline.jl")

using .CubicSpline: spline_setup, spline_eval, CubicSplineType, RealSplineType, ComplexSplineType
using .BicubicSpline: bicube_setup, bicube_eval, BicubicSplineType

export spline_setup, spline_eval, CubicSplineType, RealSplineType, ComplexSplineType
export bicube_setup, bicube_eval, BicubicSplineType

end