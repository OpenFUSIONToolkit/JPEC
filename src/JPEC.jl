module JPEC

include(joinpath(@__DIR__, "..", "deps", "build.jl"))

include("Splines/Splines.jl")
include("Vacuum/Vacuum.jl")
using .SplinesMod
using .VacuumMod


end