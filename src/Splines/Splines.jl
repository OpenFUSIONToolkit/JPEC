module SplinesMod

include("CubicSpline.jl")
include("BicubicSpline.jl")
include("FourierSpline.jl")

using .CubicSpline: spline_setup, spline_eval, spline_eval!, spline_integrate!,
 CubicSplineType, RealSplineType, ComplexSplineType
using .BicubicSpline: bicube_setup, bicube_eval, bicube_eval!, BicubicSplineType
using .FourierSpline: fspline_setup, fspline_eval, fspline_eval!, FourierSplineType

export spline_setup, spline_eval, spline_eval!, CubicSplineType, RealSplineType, ComplexSplineType
export bicube_setup, bicube_eval, bicube_eval!, BicubicSplineType
export fspline_setup, fspline_eval, fspline_eval!, FourierSplineType

end