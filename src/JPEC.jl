module JPEC

include(joinpath(@__DIR__, "..", "deps", "build.jl"))

include("Splines/Splines.jl")
using .SplinesMod


end